/*
 * Read Pileup loader - short sequencing reads stack into rows beneath a fixed
 * reference sequence (IGV-style). Matching bases are a neutral grey; mismatches
 * flare in colour. A few columns are "variant" sites where most reads agree on an
 * alt base, drawing vertical coloured stripes (an SNV) amongst the speckle of
 * sequencing noise. Reads fade in and accumulate until coverage fills, then the
 * pileup resets and rebuilds.
 */
(function () {
    "use strict";

    VGLoaders.register("pileup", {
        label: "Read Pileup",
        start: function (container, opts) {
            const BASES = "GATC", FS = 15;
            // Distinction is by hue (not brightness) so it reads on dark, light, and clear/white.
            const PAL = {
                dark:  { bg: "#07090c", ref: [110, 196, 192], match: [120, 134, 148], mismatch: [235, 120, 90] },
                light: { bg: "#eef1f5", ref: [20, 120, 116], match: [120, 132, 144], mismatch: [205, 70, 55] }
            };
            const P = PAL[opts && opts.theme === "light" ? "light" : "dark"];
            const clearBg = !!(opts && opts.clearBackground);
            const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
            const rnd = n => Math.random() * n | 0;

            const MISMATCH = 0;        // background mismatch rate (sequencing noise); 0 = only variant columns colour
            const FADE = 0.18;         // read fade-in (seconds)
            const SPAWN = 0.07;        // seconds between reads
            const GAP = 1;             // min columns between reads in a row
            const READ_LEN = 50;       // reads are a consistent length
            const HOLD = 0.9;          // seconds to hold a full pileup before resetting
            const FULL_STREAK = 25;    // failed placements in a row => pileup is full

            const cv = document.createElement("canvas");
            cv.style.width = "100%"; cv.style.height = "100%"; cv.style.display = "block";
            container.appendChild(cv);
            const ctx = cv.getContext("2d");

            let cols, rows, charW, charH, padX, padY;
            let ref, variantAlt, reads, rowEnd, spawnT, fullT, failStreak, raf, lastT;

            const rgba = (c, a) => `rgba(${c[0]},${c[1]},${c[2]},${a})`;

            function reset() {
                ref = new Uint8Array(cols);
                for (let i = 0; i < cols; i++) ref[i] = rnd(4);
                // a couple of "variant" columns with an agreed alt base
                variantAlt = new Int8Array(cols).fill(-1);
                const nVariants = Math.min(2, cols);
                for (let v = 0; v < nVariants; v++) {
                    const c = rnd(cols);
                    variantAlt[c] = (ref[c] + 1 + rnd(3)) % 4;  // guaranteed != ref
                }
                reads = [];
                rowEnd = new Array(rows - 1).fill(-GAP - 1);  // row 0 is the reference
                spawnT = 0; fullT = 0; failStreak = 0;
            }

            function init() {
                const dpr = Math.min(2, devicePixelRatio || 1);
                const w = cv.clientWidth || container.clientWidth;
                const h = cv.clientHeight || container.clientHeight;
                cv.width = w * dpr; cv.height = h * dpr;
                ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
                ctx.font = `${FS}px ${getComputedStyle(cv).fontFamily}`;
                ctx.textBaseline = "top";
                charW = ctx.measureText("G").width * 1.06; charH = FS * 1.18;
                padX = 16; padY = 14;
                cols = Math.max(8, Math.floor((w - padX * 2) / charW));
                rows = Math.max(4, Math.floor((h - padY * 2) / charH));
                reset();
            }

            // Greedy IGV-style packing: place a read in the topmost row free up to its start.
            function spawnRead() {
                const len = Math.min(READ_LEN, cols);
                const start = rnd(Math.max(1, cols - len));
                let row = -1;
                for (let r = 0; r < rowEnd.length; r++) {
                    if (rowEnd[r] < start) { row = r; break; }
                }
                if (row === -1) return false;
                rowEnd[row] = start + len + GAP;
                const mism = {};
                for (let i = 0; i < len; i++) {
                    const col = start + i;
                    if (variantAlt[col] >= 0) {
                        if (Math.random() < 0.9) mism[i] = variantAlt[col];
                    } else if (Math.random() < MISMATCH) {
                        mism[i] = (ref[col] + 1 + rnd(3)) % 4;
                    }
                }
                reads.push({ row: row + 1, start: start, len: len, mism: mism, age: 0 });
                return true;
            }

            function frame(now) {
                raf = requestAnimationFrame(frame);
                let dt = (now - (lastT || now)) / 1000; lastT = now; if (dt > 0.05) dt = 0.05;

                if (!reduce) {
                    if (fullT > 0) {
                        fullT -= dt;
                        if (fullT <= 0) reset();
                    } else {
                        spawnT -= dt;
                        while (spawnT <= 0) {
                            if (spawnRead()) {
                                failStreak = 0;
                            } else if (++failStreak > FULL_STREAK) {
                                fullT = HOLD;
                                break;
                            }
                            spawnT += SPAWN;
                        }
                    }
                    for (const rd of reads) rd.age += dt;
                }

                if (clearBg) ctx.clearRect(0, 0, cv.width, cv.height);
                else { ctx.fillStyle = P.bg; ctx.fillRect(0, 0, cv.width, cv.height); }

                // reference row
                ctx.fillStyle = rgba(P.ref, 1);
                for (let x = 0; x < cols; x++) ctx.fillText(BASES[ref[x]], padX + x * charW, padY);

                // stacked reads
                for (const rd of reads) {
                    const a = Math.min(1, rd.age / FADE);
                    const y = padY + rd.row * charH;
                    for (let i = 0; i < rd.len; i++) {
                        const mm = rd.mism[i];
                        const b = mm === undefined ? ref[rd.start + i] : mm;
                        ctx.fillStyle = rgba(mm === undefined ? P.match : P.mismatch, a);
                        ctx.fillText(BASES[b], padX + (rd.start + i) * charW, y);
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
