/*
 * Read Pulse loader - a single frozen row of GATC bases with a "reading" pulse
 * that scans L->R, flipping + glowing cells at the leading edge then locking
 * them as it trails past. Adapted from claude/loading_animations/pulse-loader.html.
 */
(function () {
    "use strict";

    VGLoaders.register("read-pulse", {
        label: "Read Pulse",
        start: function (container, opts) {
            const BASES = "GATC", FS = 15;
            // resting -> lit base colours, per theme
            const PAL = {
                dark:  { bg: "#07090c", dim: [108, 122, 136], brt: [214, 238, 234] },
                light: { bg: "#eef1f5", dim: [176, 186, 196], brt: [28, 44, 58] }
            };
            const P = PAL[opts && opts.theme === "light" ? "light" : "dark"];
            const PEAK = 0.9;          // change chance at the leading edge
            const SPEED = 23;          // cells/sec the front advances
            const WIN_FRAC = 0.5;      // pulse width as a fraction of the row (~half)
            const TAU_FRAC = 0.2;      // decay length as a fraction of the window
            const SLOW = 1;            // global playback rate (1 = original)
            const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
            const rnd = n => Math.random() * n | 0;

            const cv = document.createElement("canvas");
            cv.style.width = "100%";
            cv.style.height = "100%";
            cv.style.display = "block";
            container.appendChild(cv);
            const ctx = cv.getContext("2d");

            let cols, grid, charW, padX, rowY, WIN, TAU, head, raf, lastT;

            const COL = [];
            for (let i = 0; i <= 16; i++) {
                const v = i / 16;
                COL[i] = `rgb(${P.dim.map((d, k) => Math.round(d + (P.brt[k] - d) * v)).join(",")})`;
            }
            const col = v => COL[v <= 0 ? 0 : v >= 1 ? 16 : Math.round(v * 16)];

            function init() {
                const dpr = Math.min(2, devicePixelRatio || 1);
                const w = cv.clientWidth || container.clientWidth;
                const h = cv.clientHeight || container.clientHeight;
                cv.width = w * dpr; cv.height = h * dpr;
                ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
                ctx.font = `${FS}px ${getComputedStyle(cv).fontFamily}`;
                ctx.textBaseline = "middle";
                charW = ctx.measureText("G").width * 1.06; padX = 16;
                cols = Math.max(6, Math.floor((w - padX * 2) / charW));
                rowY = h / 2;
                WIN = Math.max(2, Math.round(cols * WIN_FRAC)); TAU = Math.max(1, WIN * TAU_FRAC);
                grid = new Uint8Array(cols); for (let i = 0; i < cols; i++) grid[i] = rnd(4);
                head = 0;
            }

            function frame(now) {
                raf = requestAnimationFrame(frame);
                let dt = (now - (lastT || now)) / 1000; lastT = now; if (dt > 0.05) dt = 0.05;
                dt *= SLOW;
                const N = cols;
                if (!reduce) {
                    head = (head + SPEED * dt) % N; const k = dt * 60;
                    const start = Math.ceil(head - WIN);
                    for (let idx = start; idx <= Math.floor(head); idx++) {
                        const i = ((idx % N) + N) % N, d = head - idx;
                        if (d >= 0 && d <= WIN && Math.random() < PEAK * Math.exp(-d / TAU) * k) grid[i] = rnd(4);
                    }
                }
                ctx.fillStyle = P.bg; ctx.fillRect(0, 0, cv.width, cv.height);
                ctx.fillStyle = col(0);
                for (let x = 0; x < cols; x++) ctx.fillText(BASES[grid[x]], padX + x * charW, rowY);
                if (!reduce) {
                    const start = Math.ceil(head - WIN);
                    for (let idx = start; idx <= Math.floor(head); idx++) {
                        const i = ((idx % N) + N) % N, d = head - idx; if (d < 0 || d > WIN) continue;
                        const v = Math.exp(-d / TAU); ctx.fillStyle = col(v);
                        ctx.fillText(BASES[grid[i]], padX + i * charW, rowY);
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
        }
    });
})();
