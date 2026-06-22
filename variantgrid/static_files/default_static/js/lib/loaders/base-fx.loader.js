/*
 * Base FX loaders - ambient brightness patterns over an anonymous GATC field.
 * Three modes: matrix (falling drops), ripple (expanding rings), wave.
 * Adapted from claude/loading_animations/base-fx-loader.html. Each mode is
 * registered as its own selectable loader via a shared factory.
 */
(function () {
    "use strict";

    function makeBaseFx(MODE) {
        return function (container, opts) {
            const BASES = "GATC", FS = 15;
            const PAL = {
                dark:  { bg: "#07090c", dim: [108, 122, 136], brt: [158, 182, 186] },
                light: { bg: "#eef1f5", dim: [182, 190, 198], brt: [70, 96, 112] }
            };
            const P = PAL[opts && opts.theme === "light" ? "light" : "dark"];
            const clearBg = !!(opts && opts.clearBackground);   // transparent canvas (shows page bg)
            const SLOW = 1;     // global playback rate (1 = original)
            const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
            const rnd = n => Math.random() * n | 0;

            const cv = document.createElement("canvas");
            cv.style.width = "100%";
            cv.style.height = "100%";
            cv.style.display = "block";
            container.appendChild(cv);
            const ctx = cv.getContext("2d");

            let cols, rows, grid, cool, charW, charH, drops = [], spawnT = 0, t0 = 0, raf, lastT, pad;

            const COL = [];
            for (let i = 0; i <= 16; i++) {
                const v = i / 16;
                COL[i] = `rgb(${P.dim.map((d, k) => Math.round(d + (P.brt[k] - d) * v)).join(",")})`;
            }
            const col = v => COL[Math.max(0, Math.min(16, Math.round(v * 16)))];

            function init() {
                const dpr = Math.min(2, devicePixelRatio || 1);
                const w = cv.clientWidth || container.clientWidth;
                const h = cv.clientHeight || container.clientHeight;
                cv.width = w * dpr; cv.height = h * dpr;
                ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
                ctx.font = `${FS}px ${getComputedStyle(cv).fontFamily}`;
                ctx.textBaseline = "top";
                charW = ctx.measureText("G").width * 1.06; charH = FS * 1.02;
                pad = 16;
                cols = Math.max(6, Math.floor((w - pad * 2) / charW));
                rows = Math.max(4, Math.floor((h - pad * 2) / charH));
                grid = new Uint8Array(cols * rows); for (let i = 0; i < grid.length; i++) grid[i] = rnd(4);
                cool = new Float32Array(cols * rows);   // ripple: per-cell churn cooldown (seconds)
                drops = [];
            }

            function intensity(I, time, dt) {
                if (MODE === "matrix") {
                    spawnT -= dt;
                    if (spawnT <= 0 && drops.length < cols * 0.55) {
                        drops.push({ x: rnd(cols), y: -rnd(6), v: 6 + Math.random() * 10, t: 8 + rnd(11) });
                        spawnT = 0.09 + Math.random() * 0.14;
                    }
                    for (const d of drops) d.y += d.v * dt;
                    drops = drops.filter(d => d.y - d.t < rows);
                    for (const d of drops) {
                        const head = d.y;
                        for (let y = 0; y < rows; y++) {
                            if (y <= head) {
                                const dd = head - y;
                                if (dd < d.t) { let v = 1 - dd / d.t; if (y === (head | 0)) v = 1; const i = y * cols + d.x; if (v > I[i]) I[i] = v; }
                            }
                        }
                    }
                } else if (MODE === "ripple") {
                    spawnT -= dt;
                    if (spawnT <= 0) {
                        drops.push({ cx: 4 + Math.random() * (cols - 8), cy: 3 + Math.random() * (rows - 6), age: 0 });
                        spawnT = 0.97 + Math.random() * 0.75;
                    }
                    // Ring grows (in cells) to the full grid diagonal over its life so it
                    // always reaches the edge, staying bright until it gets there then fading.
                    const LIFE = 7.5, maxR = Math.hypot(cols, rows);
                    for (const d of drops) d.age += dt;
                    drops = drops.filter(d => d.age < LIFE);
                    for (const d of drops) {
                        const p = d.age / LIFE, r = p * maxR;        // r in cells
                        const dec = p < 0.8 ? 1 : Math.max(0, (1 - p) / 0.2);
                        for (let y = 0; y < rows; y++) for (let x = 0; x < cols; x++) {
                            const dist = Math.hypot((x - d.cx) * charW, (y - d.cy) * charH) / charW;
                            const v = Math.exp(-(((dist - r) / 1.7) ** 2)) * dec;
                            if (v > 0.02) { const i = y * cols + x; I[i] = Math.min(1, I[i] + v); }
                        }
                    }
                } else { // wave
                    const wl = cols / 1.6, ph0 = time * 1.9;
                    for (let y = 0; y < rows; y++) for (let x = 0; x < cols; x++) {
                        const ph = (x - y * 0.22) * (2 * Math.PI / wl) - ph0; const s = (Math.sin(ph) + 1) / 2;
                        I[y * cols + x] = s * s;
                    }
                }
            }

            function frame(now) {
                raf = requestAnimationFrame(frame);
                if (!t0) t0 = now; const time = (now - t0) / 1000 * SLOW;
                let dt = (now - (lastT || now)) / 1000; lastT = now; if (dt > 0.05) dt = 0.05;
                dt *= SLOW;

                const I = new Float32Array(cols * rows);
                if (!reduce) intensity(I, time, dt);

                if (clearBg) ctx.clearRect(0, 0, cv.width, cv.height);
                else { ctx.fillStyle = P.bg; ctx.fillRect(0, 0, cv.width, cv.height); }
                for (let y = 0; y < rows; y++) {
                    const py = pad + y * charH;
                    for (let x = 0; x < cols; x++) {
                        const i = y * cols + x, v = I[i];
                        if (MODE === "ripple") {
                            // The ripple "wakes" a cell as it passes; it keeps rerolling for a
                            // short random cooldown afterwards, then freezes again.
                            if (v > 0.15) cool[i] = 0.3 + Math.random() * 0.9;
                            if (cool[i] > 0) {
                                cool[i] -= dt;
                                if (Math.random() < 0.04 * (dt * 60)) grid[i] = rnd(4);
                            }
                        } else if (Math.random() < (0.012 + 0.45 * v) * (dt * 60)) {
                            grid[i] = rnd(4);   // matrix/wave: churn scaled by brightness
                        }
                        ctx.fillStyle = col(v);
                        ctx.fillText(BASES[grid[i]], pad + x * charW, py);
                    }
                }
            }

            let resizeTimer;
            const onResize = () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(init, 150); };
            addEventListener("resize", onResize);

            init();
            raf = requestAnimationFrame(frame);

            return function stop() {
                cancelAnimationFrame(raf);
                removeEventListener("resize", onResize);
                clearTimeout(resizeTimer);
                cv.remove();
            };
        };
    }

    VGLoaders.register("base-fx-matrix", { label: "Base FX - Matrix", start: makeBaseFx("matrix") });
    VGLoaders.register("base-fx-ripple", { label: "Base FX - Ripple", start: makeBaseFx("ripple") });
    VGLoaders.register("base-fx-wave", { label: "Base FX - Wave", start: makeBaseFx("wave") });
})();
