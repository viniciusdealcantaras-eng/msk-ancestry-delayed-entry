#!/usr/bin/env Rscript
# Figuras do projeto Ancestralidade MSK (base R + survival). Uso: Rscript figuras.R
suppressPackageStartupMessages(library(survival))
od <- "out/figs"; dir.create(od, showWarnings = FALSE)
FONT <- "Liberation Sans"; CUR <- FONT
save_fig <- function(name, w, h, fun) {
  CUR <<- FONT
  tiff(file.path(od, paste0(name, ".tif")), width = w, height = h, units = "in", res = 600,
       compression = "lzw", type = "cairo", family = FONT); fun(); dev.off()
  png(file.path(od, paste0(name, ".png")), width = w, height = h, units = "in", res = 300,
      type = "cairo", family = FONT); fun(); dev.off()
  CUR <<- "Helvetica"
  pdf(file.path(od, paste0(name, ".pdf")), width = w, height = h, family = "Helvetica",
      useDingbats = FALSE); fun(); dev.off()
}
GL <- c("nonASJ-EUR" = "European (non-Ashkenazi)", "ASJ-EUR" = "European (Ashkenazi)",
        "ADMIX" = "Admixed/other", "AFR" = "African", "EAS" = "East Asian",
        "SAS" = "South Asian", "NAM" = "Native American")
ORD <- names(GL)
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# ---------------- Figure 1: composicao por tipo de cancer ----------------
cp <- read.csv("out/fig_comp.csv", check.names = FALSE, row.names = 1)
cp <- cp[, ORD]; cp <- cp[order(rowSums(cp)), ]
pr <- cp / rowSums(cp)
shades <- c("grey92", "grey78", "grey62", "black", "grey40", "grey25", "white")
dens <- c(NA, NA, NA, NA, NA, NA, 30)
f1 <- function() {
  par(mar = c(7.2, 13.5, 0.8, 4.2), family = CUR, lwd = 0.7, xpd = NA)
  b <- barplot(t(as.matrix(pr)), horiz = TRUE, col = shades, border = "grey20", space = 0.35,
               axes = FALSE, names.arg = rep("", nrow(pr)), lwd = 0.5)
  axis(1, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20)), cex.axis = 0.8,
       mgp = c(3, 0.45, 0), tcl = -0.3, xpd = FALSE)
  mtext("Patients (%)", side = 1, line = 1.6, cex = 0.8)
  text(-0.02, b, rownames(pr), adj = 1, cex = 0.78)
  text(1.02, b, format(rowSums(cp), big.mark = ","), adj = 0, cex = 0.72)
  text(1.02, max(b) + 0.95, "n", adj = 0, cex = 0.72, font = 2)
  lx <- c(-0.27, -0.27, 0.2, 0.2, 0.48, 0.48, 0.76); ly <- c(-4.1, -5.0, -4.1, -5.0, -4.1, -5.0, -4.1)
  for (i in 1:7) { rect(lx[i], ly[i] - 0.22, lx[i] + 0.03, ly[i] + 0.22, col = shades[i], border = "grey20", lwd = 0.5)
                   text(lx[i] + 0.045, ly[i], GL[i], adj = 0, cex = 0.7) }
}
save_fig("Figure1", 7.2, 5.6, f1)

# ---------------- Figure 2: intervalo diagnostico-sequenciamento ----------------
iv <- read.csv("out/interval.csv")
qr <- read.csv("out/fig_qreg.csv")
G6 <- c("nonASJ-EUR", "ASJ-EUR", "ADMIX", "AFR", "EAS", "SAS")
lt <- c(1, 2, 3, 1, 4, 5); lc <- c("grey55", "black", "black", "black", "black", "black")
lw <- c(2.2, 1, 1, 1.8, 1, 1)
f2 <- function() {
  layout(matrix(1:2, 1), widths = c(1.15, 1))
  par(mar = c(3.6, 4, 1.8, 0.8), family = CUR, lwd = 0.8)
  plot(NA, xlim = log10(c(0.1, 300)), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_along(G6)) {
    x <- sort(pmax(iv$dx_to_seq_months[iv$anc == G6[i]], 0.1))
    lines(log10(c(0.1, x)), c(0, seq_along(x) / length(x)), type = "s", lty = lt[i], col = lc[i], lwd = lw[i])
  }
  tk <- c(0.1, 1, 3, 6, 12, 24, 60, 120, 240)
  axis(1, at = log10(tk), labels = tk, cex.axis = 0.78, mgp = c(3, 0.45, 0), tcl = -0.3)
  axis(2, at = seq(0, 1, 0.2), labels = fmt(seq(0, 1, 0.2), 1), las = 1, cex.axis = 0.8,
       mgp = c(3, 0.55, 0), tcl = -0.3)
  box(bty = "l")
  mtext("Months from diagnosis to sequencing (log scale)", side = 1, line = 1.8, cex = 0.78)
  mtext("Cumulative proportion sequenced", side = 2, line = 2.6, cex = 0.78)
  mtext("A", side = 3, line = 0.4, adj = -0.18, font = 2)
  legend("bottomright", legend = GL[G6], lty = lt, col = lc, lwd = lw, bty = "n", cex = 0.68, seg.len = 2.6)
  # painel B
  par(mar = c(3.6, 7.6, 1.8, 1.6))
  q <- qr[match(rev(G6[-1]), qr$grp), ]; y <- seq_len(nrow(q))
  xl <- c(-3, 5.5)
  plot(NA, xlim = xl, ylim = c(0.5, nrow(q) + 0.5), axes = FALSE, xlab = "", ylab = "")
  abline(v = 0, lty = 3, col = "grey40")
  segments(q$clo, y + 0.15, q$chi, y + 0.15, col = "grey50"); points(q$crude, y + 0.15, pch = 1, cex = 0.85)
  segments(q$alo, y - 0.15, q$ahi, y - 0.15); points(q$adj, y - 0.15, pch = 16, cex = 0.85)
  axis(1, at = seq(-2, 5, 1), cex.axis = 0.78, mgp = c(3, 0.45, 0), tcl = -0.3)
  axis(2, at = y, labels = GL[q$grp], las = 1, cex.axis = 0.75, tick = FALSE, mgp = c(3, 0.3, 0))
  box(bty = "l")
  mtext("Difference in median (months)", side = 1, line = 1.8, cex = 0.78)
  mtext("B", side = 3, line = 0.4, adj = -0.75, font = 2)
  legend("bottomright", legend = c("Unadjusted", "Adjusted"), pch = c(1, 16), bty = "n", cex = 0.68)
}
save_fig("Figure2", 7.2, 3.6, f2)

# ---------------- Figure 3: Kaplan-Meier ingenuo x corrigido ----------------
s <- read.csv("out/surv_dx.csv")
G3 <- c("nonASJ-EUR", "AFR", "EAS"); lt3 <- c(1, 1, 2); lc3 <- c("grey55", "black", "black"); lw3 <- c(2.2, 1.6, 1.2)
times <- seq(0, 120, 24)
f3 <- function() {
  layout(matrix(1:2, 1))
  for (p in 1:2) {
    par(mar = c(8.4, 8.6, 1.8, 1.2), family = CUR, lwd = 0.8)
    plot(NA, xlim = c(0, 120), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    cnt <- list()
    for (i in seq_along(G3)) {
      x <- s[s$anc == G3[i], ]
      f <- if (p == 1) survfit(Surv(T, E) ~ 1, data = x) else survfit(Surv(L, T, E) ~ 1, data = x)
      lines(f, conf.int = FALSE, lty = lt3[i], col = lc3[i], lwd = lw3[i], mark.time = FALSE)
      cnt[[i]] <- sapply(times, function(t) if (p == 1) sum(x$T >= t) else sum(x$L <= t & x$T >= t))
    }
    axis(1, at = times, cex.axis = 0.8, mgp = c(3, 0.45, 0), tcl = -0.3)
    axis(2, at = seq(0, 1, 0.2), labels = fmt(seq(0, 1, 0.2), 1), las = 1, cex.axis = 0.8, mgp = c(3, 0.55, 0), tcl = -0.3)
    box(bty = "l")
    mtext("Months from diagnosis", side = 1, line = 1.7, cex = 0.8)
    mtext("Overall survival", side = 2, line = 2.4, cex = 0.8)
    mtext(if (p == 1) "A  Without delayed entry" else "B  With delayed entry", side = 3, line = 0.4, adj = -0.32, font = 2, cex = 0.85)
    mtext("Number at risk", side = 1, line = 3.2, at = -15, adj = 1, cex = 0.66, font = 2)
    for (i in seq_along(G3)) {
      mtext(c("European", "African", "East Asian")[i], side = 1, line = 3.2 + i * 0.95, at = -15, adj = 1, cex = 0.66)
      mtext(format(cnt[[i]], big.mark = ","), side = 1, line = 3.2 + i * 0.95, at = times, cex = 0.66)
    }
    if (p == 1) legend("bottomleft", legend = c("European (non-Ashk.)", "African", "East Asian"), lty = lt3, col = lc3, lwd = lw3, bty = "n", cex = 0.66, seg.len = 2.2)
  }
}
save_fig("Figure3", 7.2, 4.0, f3)

# ---------------- Figure 4: forest de HR ----------------
fo <- read.csv("out/fig_forest.csv")
fo <- fo[fo$grp != "NAM", ]
mods <- unique(fo$model); pch4 <- c(1, 16, 15, 2, 17, 18)
f4 <- function() {
  par(mar = c(3.4, 12.5, 0.8, 9.5), family = CUR, lwd = 0.8, xpd = NA)
  G5 <- c("ASJ-EUR", "ADMIX", "AFR", "EAS", "SAS"); nm <- length(mods)
  yy <- c(); lab <- c(); y0 <- 0
  plot(NA, xlim = log(c(0.6, 1.7)), ylim = c(0.3, length(G5) * (nm + 1.2)), axes = FALSE, xlab = "", ylab = "")
  segments(0, 0.3, 0, length(G5) * (nm + 1.2) - 0.5, lty = 3, col = "grey40")
  ytop <- length(G5) * (nm + 1.2)
  for (j in seq_along(G5)) {
    base <- ytop - (j - 1) * (nm + 1.2)
    text(log(0.6) - 0.02, base - 0.4, GL[G5[j]], adj = 1, font = 2, cex = 0.74)
    for (k in seq_along(mods)) {
      r <- fo[fo$grp == G5[j] & fo$model == mods[k], ]; y <- base - 0.4 - k
      segments(log(r$lo), y, log(r$hi), y)
      points(log(r$hr), y, pch = pch4[k], cex = 0.8, bg = "white")
      text(log(0.6) - 0.02, y, mods[k], adj = 1, cex = 0.64)
      text(log(1.7) + 0.04, y, paste0(fmt(r$hr), " (", fmt(r$lo), "-", fmt(r$hi), ")"), adj = 0, cex = 0.64)
    }
  }
  text(log(1.7) + 0.04, ytop + 0.2, "HR (95% CI)", adj = 0, cex = 0.68, font = 2)
  tk <- c(0.6, 0.8, 1, 1.25, 1.5)
  axis(1, at = log(tk), labels = tk, cex.axis = 0.78, mgp = c(3, 0.45, 0), tcl = -0.3)
  mtext("Hazard ratio vs European (non-Ashkenazi), log scale", side = 1, line = 1.7, cex = 0.78)
}
save_fig("Figure4", 7.2, 7.4, f4)

# ---------------- Figure 5: por tipo de cancer ----------------
pc <- read.csv("out/fig_percancer.csv")
fS1 <- function() {
  par(mar = c(3.4, 13.5, 0.8, 1), family = CUR, lwd = 0.8, xpd = NA)
  keys <- unique(pc[, c("ct", "grp")]); keys <- keys[order(keys$ct, keys$grp), ]
  n <- nrow(keys); plot(NA, xlim = log(c(0.6, 2.4)), ylim = c(0.4, n + 1.3), axes = FALSE, xlab = "", ylab = "")
  segments(0, 0.4, 0, n + 0.4, lty = 3, col = "grey40")
  for (i in seq_len(n)) {
    y <- n + 1 - i
    for (m in c("naive", "corr")) {
      r <- pc[pc$ct == keys$ct[i] & pc$grp == keys$grp[i] & pc$model == m, ]; dy <- if (m == "naive") 0.14 else -0.14
      segments(log(r$lo), y + dy, log(r$hi), y + dy, col = if (m == "naive") "grey50" else "black")
      points(log(r$hr), y + dy, pch = if (m == "naive") 1 else 16, cex = 0.8)
    }
    text(log(0.6) - 0.03, y, paste0(keys$ct[i], ", ", ifelse(keys$grp[i] == "AFR", "African", "East Asian"),
         " (n = ", r$n_grp, ")"), adj = 1, cex = 0.68)
  }
  tk <- c(0.6, 0.8, 1, 1.5, 2); axis(1, at = log(tk), labels = tk, cex.axis = 0.78, mgp = c(3, 0.45, 0), tcl = -0.3)
  mtext("Hazard ratio vs European (non-Ashkenazi), log scale", side = 1, line = 1.7, cex = 0.78)
  legend(x = log(0.6), y = n + 1.45, legend = c("Without delayed entry", "With delayed entry"), pch = c(1, 16),
         col = c("grey50", "black"), bty = "n", cex = 0.68, horiz = TRUE)
}
save_fig("Figure5", 7.2, 4.6, fS1)
cat("ok\n")
