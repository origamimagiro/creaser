import { M } from "./flatfolder/math.js";
import { NOTE } from "./flatfolder/note.js";
import { SVG } from "./flatfolder/svg.js";
import { IO } from "./flatfolder/io.js";
import { X } from "./flatfolder/conversion.js";
import { SOLVER } from "./flatfolder/solver.js";
import { CON } from "./flatfolder/constraints.js";

window.onload = () => { MAIN.startup(); };  // entry point

const MAIN = {
    color: {
        background: "lightgray",
        paper: "white",
        B: "black",
        C: "green",
        R: "red",
    },
    size: {
        b: 50,
        s: SVG.SCALE,
    },
    startup: () => {
        CON.build();
        NOTE.clear_log();
        NOTE.start("*** Starting Creaser ***");
        NOTE.time("Initializing interface");
        const [b, s] = [MAIN.size.b, MAIN.size.s];
        const main = document.getElementById("main");
        for (const [k, v] of Object.entries({
            xmlns: SVG.NS,
            style: `background: ${MAIN.color.background}`,
            viewBox: [0, 0, 2*s, s].join(" "),
        })) {
            main.setAttribute(k, v);
        }
        for (const [i, id] of ["input", "output"].entries()) {
            const svg = document.getElementById(id);
            for (const [k, v] of Object.entries({
                xmlns: SVG.NS,
                height: s,
                width: s,
                x: i*s,
                y: 0,
                viewBox: [-b, -b, s + 2*b, s + 2*b].join(" "),
            })) {
                svg.setAttribute(k, v);
            }
        }
        document.getElementById("import").onchange = (e) => {
            if (e.target.files.length > 0) {
                const file_reader = new FileReader();
                file_reader.onload = MAIN.process_file;
                file_reader.readAsText(e.target.files[0]);
            }
        };
    },
    process_file: (e) => {
        NOTE.clear_log();
        NOTE.start("*** Starting File Import ***");
        const doc = e.target.result;
        const file_name = document.getElementById("import").value;
        const parts = file_name.split(".");
        const type = parts[parts.length - 1].toLowerCase();
        if (type != "fold") {
            console.log(`Found file with extension ${type}, FOLD format required`);
            return;
        }
        NOTE.time(`Importing from file ${file_name}`);
        const ex = JSON.parse(doc);
        const properties = ["vertices_coords", "edges_vertices", "edges_assignment"];
        const [V_org, EVi, EAi] = properties.map(property => {
            const val = ex[property];
            if (val == undefined) {
                NOTE.time(`FOLD file must contain ${property}, but not found`);
                return undefined;
            }
            return val;
        });
        for (const property of [V_org, EVi]) {
            if (property == undefined) { return; }
        }
        const Vi = M.normalize_points(V_org);
        const L = EVi.map((P) => M.expand(P, Vi));
        const eps = M.min_line_length(L) / M.EPS;
        const [C, VC] = MAIN.V_2_C_VC(Vi, eps);
        const Xo = [];
        const Do = [];
        for (let i = 0; i < EVi.length; ++i) {
            const [v1, v2] = EVi[i].map(v => VC[v].map(c => C[c]));
            const [[x, y], d] = MAIN.coords_2_line(v1, v2);
            Xo.push([x, i]);
            Do.push([d, i]);
        }
        Xo.sort(([x1, i1], [x2, i2]) => x1 - x2);
        Do.sort(([d1, i1], [d2, i2]) => d1 - d2);
        const EL = EVi.map(() => [undefined, undefined]);
        const X = [Xo[0][0]];
        EL[Xo[0][1]][0] = 0;
        for (let i = 1; i < Xo.length; ++i) {
            const [x1, i1] = Xo[i - 1];
            const [x2, i2] = Xo[i];
            if (Math.abs(x2 - x1) > 10000*M.FLOAT_EPS) {
                X.push(x2);
            }
            EL[i2][0] = X.length - 1;
        }
        const D = [Do[0][0]];
        EL[Do[0][1]][1] = 0;
        for (let i = 1; i < Do.length; ++i) {
            const [d1, i1] = Do[i - 1];
            const [d2, i2] = Do[i];
            if (Math.abs(d2 - d1) > 10000*M.FLOAT_EPS) {
                D.push(d2);
            }
            EL[i2][1] = D.length - 1;
        }
        const EPS = 10000*M.FLOAT_EPS;
        const ES = EL.map(([xi, di]) => (
            ((Math.abs(X[xi]) < EPS) ||
            (Math.abs(X[xi] - 1) < EPS)) &&
            ((Math.abs(D[di]) < EPS) ||
            (Math.abs(D[di] - 1) < EPS))
        ));
        const Y = X.map(x => (1 - x*x)**0.5);
        const target = {C, VC, X, Y, D, EL, ES, EV: EVi, EA: EAi};
        MAIN.update(target, SVG.clear("input"));
        const FOLD = {
            C: [0, 1],
            VC: [[0, 0], [0, 1], [1, 0], [1, 1]],
            EV: [[0, 1], [0, 2], [1, 3], [2, 3]],
            EA: ["B", "B", "B", "B"],
            ES: [true, true, true, true],
        };
        const FS = [FOLD];
        MAIN.update(FOLD, SVG.clear("output"));
        const prev = document.getElementById("prev");
        let i = 0;
        prev.onclick = () => {
            if (i > 0) {
                --i;
                MAIN.update(FS[i], SVG.clear("output"));
                MAIN.update(target, SVG.clear("input"));
            }
        };
        const next = document.getElementById("next");
        next.onclick = () => {
            ++i;
            if (i < FS.length) {
                MAIN.update(FS[i], SVG.clear("output"));
                MAIN.update(target, SVG.clear("input"));
            } else if (i == FS.length) {
                console.log("looking for next line...");
                const FOLD2 = MAIN.get_next(target, FS[i - 1], eps);
                if (FOLD2 == undefined) {
                    console.log("no line found!");
                    --i;
                } else {
                    console.log("found a new line");
                    FS.push(FOLD2);
                    console.log(FOLD2);
                    MAIN.update(FOLD2, SVG.clear("output"));
                    MAIN.update(target, SVG.clear("input"));
                }
            }
        };
    },
    V_2_C_VC: (V, eps) => {
        const Ci = [];
        for (let i = 0; i < V.length; ++i) {
            for (const j of [0, 1]) {
                Ci.push([V[i][j], i, j]);
            }
        }
        Ci.sort(([a, ai, aj], [b, bi, bj]) => a - b);
        const C = [];
        const VC = V.map(() => [undefined, undefined]);
        C.push(Ci[0][0]);
        VC[Ci[0][1]][Ci[0][2]] = 0;
        for (let i = 1; i < Ci.length; ++i) {
            const [c1, i1, j1] = Ci[i - 1];
            const [c2, i2, j2] = Ci[i];
            if (c2 - c1 > eps) {
                C.push(c2);
            }
            VC[i2][j2] = C.length - 1;
        }
        return [C, VC];
    },
    get_idx: (C, c, eps = 10000*M.FLOAT_EPS) => {
        let [l, r] = [0, C.length];
        while (r - l > 1) {
            const m = (r + l) >> 1;
            const diff = C[m] - c;
            if (Math.abs(diff) < eps) {
                return m;
            }
            if (diff > 0) {
                r = m;
            } else {
                l = m;
            }
        }
        if (l == C.length) {
            return undefined;
        }
        return (Math.abs(C[l] - c) < eps) ? l : undefined;
    },
    get_next: (target, FOLD, eps) => {
        const Ci = target.C;
        const VCi = target.VC;
        const EVi = target.EV;
        const {EL, ES, X, Y, D} = target;
        const {C, VC, EV, EA} = FOLD;
        const V = VC.map(c => M.expand(c, C));
        const Vi = VCi.map(c => M.expand(c, Ci));
        const EL_map = new Map();
        for (let i = 0; i < EL.length; ++i) {
            const [xi, di] = EL[i];
            const s = ES[i];
            let out = EL_map.get(di);
            if (out == undefined) {
                out = new Map();
                EL_map.set(di, out);
            }
            out.set(xi, s);
        }
        let found = undefined;
        const L = MAIN.get_line(V, EL_map, D, X);
        if (L == undefined) { return undefined; }
        const [xi, di, P] = L;
        for (let i = 0; i < EL.length; ++i) {
            if (ES[i]) { continue; }
            const [xj, dj] = EL[i];
            if ((xi == xj) && (di == dj)) {
                ES[i] = true;
            }
        }
        const line = [[X[xi], Y[xi]], D[di]];
        const [EV2, EA2, ES2, V2] = MAIN.EV_EA_V_line_eps_2_EV2_EA2_ES2_V2(EV, EA, V, line);
        const [C2, VC2] = MAIN.V_2_C_VC(V2, eps);
        const FOLD2 = {C: C2, VC: VC2, EV: EV2, EA: EA2, ES: ES2, P};
        return FOLD2;
    },
    update: (FOLD, svg) => {
        const {C, VC, EV, EA, ES, P} = FOLD;
        SVG.append("rect", svg, {
            x: 0,
            y: 0,
            width: MAIN.size.s,
            height: MAIN.size.s,
            fill: MAIN.color.paper,
        });
        const V = VC.map(c => M.expand(c, C));
        const lines = EV.map(l => M.expand(l, V));
        SVG.draw_segments(svg, lines, {
            id: "flat_e_boundary", stroke: MAIN.color.C,
            stroke_width: 1, filter: i => !ES[i],
        });
        SVG.draw_segments(svg, lines, {
            id: "flat_e_boundary", stroke: MAIN.color.R,
            stroke_width: 1, filter: i => ES[i],
        });
        if (P != undefined) {
            SVG.draw_points(svg, P.map(i => V[i]), {
                id: "flat_p", fill: MAIN.color.R, r: 10,
            });
        }
    },
    get_line: (P, EL_map, D, X) => {
        const check = ([[x, y], d]) => {
            const di = MAIN.get_idx(D, d);
            if (di == undefined) { return; }
            const Xi = EL_map.get(di);
            if (Xi == undefined) { return; }
            const xi = MAIN.get_idx(X, x);
            if (xi == undefined) { return; }
            const s = Xi.get(xi);
            if (s == undefined) { return; }
            if (s == true) { return; }
            return [xi, di];
        };
        const n = P.length;
        for (let i1 = 0; i1 < n; ++i1) {
            const a = P[i1];
            for (let i2 = i1 + 1; i2 < n; ++i2) {
                const b = P[i2];
                const m = M.div(M.add(a, b), 2);
                const v1 = M.sub(b, a); // line through AB
                const v2 = M.perp(v1);  // perpendicular bisector of AB
                for (let v of [v1, v2]) {
                    const u = M.unit(v);
                    const d = M.dot(m, u);
                    const L = MAIN.normalize_line([u, d]);
                    const out = check(L);
                    if (out != undefined) {
                        out.push([i1, i2]);
                        return out;
                    }
                }
            }
        }
        for (let i1 = 0; i1 < n; ++i1) {
            const a = P[i1];
            for (let i2 = i1 + 1; i2 < n; ++i2) {
                const b = P[i2];
                for (let i3 = i2 + 1; i3 < n; ++i3) {
                    const c = P[i3];
                    for (const [a1, b1, c1] of [
                        [a, b, c],
                        [c, a, b],
                        [b, c, a],
                    ]) {
                        if (Math.abs(M.area2(a1, b1, c1)) > 0.000001) {
                            {   // angle bisector of ABC
                                const ba = M.unit(M.sub(a1, b1));
                                const bc = M.unit(M.sub(c1, b1));
                                const u = M.unit(M.add(ba, bc));
                                const v = M.perp(u);
                                const d = M.dot(b1, v);
                                const L = MAIN.normalize_line([v, d]);
                                const out = check(L);
                                if (out != undefined) {
                                    out.push([i1, i2, i3]);
                                    return out;
                                }
                            }
                            {   // through A perpendicular to BC
                                const u = M.unit(M.sub(c1, b1));
                                const d = M.dot(a1, u);
                                const L = MAIN.normalize_line([u, d]);
                                const out = check(L);
                                if (out != undefined) {
                                    out.push([i1, i2, i3]);
                                    return out;
                                }
                            }
                        } else {
                            {   // through A perpendicular to BC
                                const u = M.unit(M.sub(c1, b1));
                                const d = M.dot(a1, u);
                                const L = MAIN.normalize_line([u, d]);
                                const out = check(L);
                                if (out != undefined) {
                                    out.push([i1, i2, i3]);
                                    return out;
                                }
                            }
                        }
                    }
                }
            }
        }
    },
    normalize_line: ([u, d]) => {
        let out = [u, d];
        if ((u[0] + 1) < M.FLOAT_EPS)  {
            out = [M.mul(u, -1), -d];
        }
        if (u[1] < -M.FLOAT_EPS)  {
            out = [M.mul(u, -1), -d];
        }
        return out;
    },
    coords_2_line: (a, b) => {
        const v = M.sub(b, a);
        const u = M.perp(M.unit(v));
        const d = M.dot(a, u);
        return MAIN.normalize_line([u, d]);
    },
    line_2_coords: (line) => {
        const [u, d] = line;
        const p = M.mul(u, d);
        const off = M.mul(M.perp(u), 10);
        const p1 = M.add(p, off);
        const p2 = M.sub(p, off);
        return [p1, p2];
    },
    EV_EA_V_line_eps_2_EV2_EA2_ES2_V2: (EV, EA, V, line, eps = M.FLOAT_EPS) => {
        // assumes convex faces (or line divides a face into at most two pieces
        const [u, d] = line;
        const [a, b] = MAIN.line_2_coords(line);
        const nV = V.length;
        const V2 = V.map(v => v);
        const VD = V.map(v => {
            const dv = M.dot(u, v) - d;
            return (Math.abs(dv) <= eps) ? 0 : dv;
        });
        const EVA = EV.map(([v1, v2], i) => [v1, v2, EA[i], false]);
        const P = new Set();
        for (let i = 0; i < EV.length; ++i) {
            const [v1, v2, ea, s] = EVA[i];
            const on1 = (Math.abs(VD[v1]) < 10000*eps);
            const on2 = (Math.abs(VD[v2]) < 10000*eps);
            if (on1) { P.add(v1); }
            if (on2) { P.add(v2); }
            if (on1 || on2) { continue; }
            if ((VD[v1] < 0) != (VD[v2] < 0)) {
                const x = M.intersect([V[v1], V[v2]], [a, b], eps);
                if (x == undefined) { continue; }
                const xi = V2.length;
                V2.push(x);
                P.add(xi);
                EVA[i][1] = xi;
                EVA.push([xi, v2, ea, false]);
            }
        }
        const perp = M.perp(u);
        const A = Array.from(P);
        A.sort((v1, v2) => M.dot(perp, V2[v1]) - M.dot(perp, V2[v2]));
        for (let i = 1; i < A.length; ++i) {
            EVA.push([A[i - 1], A[i], "F", true]);
        }
        for (let i = 1; i < EVA.length; ++i) {
            let [v1, v2, ea, s] = EVA[i];
            if (v2 < v1) {
                [v1, v2] = [v2, v1];
            }
            EVA[i] = [v1, v2, ea, s];
        }
        EVA.sort((va1, va2) => {
            const diff = va1[0] - va2[0];
            return (diff == 0) ? (va1[1] - va2[1]) : diff;
        });
        const EV2 = EVA.map(([v1, v2, ea, s]) => [v1, v2]);
        const EA2 = EVA.map(([v1, v2, ea, s]) => ea);
        const ES2 = EVA.map(([v1, v2, ea, s]) => s);
        return [EV2, EA2, ES2, V2];
    },
    FV_VD_2_FG: (FV, VD) => {
        const EF_map = new Map();
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                EF_map.set(M.encode([v2, v1]), i);
            }
        }
        const FG = FV.map(() => undefined);
        let g = 0;
        const dfs = (i) => {
            if (FG[i] != undefined) { return; }
            FG[i] = g;
            const F = FV[i];
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                if ((VD[v1] == 0) && (VD[v2] == 0)) { continue; }
                const fi = EF_map.get(M.encode([v1, v2]));
                if (fi != undefined) {
                    dfs(fi);
                }
            }
        };
        for (let i = 0; i < FG.length; ++i) {
            if (FG[i] != undefined) { continue; }
            dfs(i);
            ++g;
        }
        return FG;
    },
};
