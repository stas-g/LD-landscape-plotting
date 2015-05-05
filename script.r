#code for plotting landscape LD with a comparison heatmap
#requires library 'plotrix'
#v1.2 --> nfg, 08.04.2015

print('landscape_script: v1.2, 08.04.2015')

#ld.win --> for windows of a specified width n (in number of SNPs) calculate and return a vector of average LDs for each window. if SNP positions (pos) are provided then output is a matrix with first column = window LDs and second column = window coordinates (position of the middle SNP in a window, if position). if plot_ = TRUE, landscape plot of window LDs against (real) window positions is also plotted. 
# 18.03.2015 --> v1.1

ld.win <- function(n, r2mat, pos = NULL, plot_ = FALSE, ...){
		no.win <- nrow(r2mat) - n + 1
		ans <- sapply(1 : no.win, FUN = function(i){
				ind <- (1 : n) + i - 1
				R2.win <- mean(r2mat[ind, ind], na.rm = TRUE)
				if(!is.null(pos)){
					pos.win <- pos[floor(n/2) + i - 1]
					return(c(R2.win = R2.win, pos.win = pos.win))
						} else { return(R2.win) }
			})
		if(!is.null(pos)) ans <- t(ans)
		if(plot_ == TRUE & !is.null(pos)) plot(ans[, 'pos.win'], ans[, 'R2.win'], ...)
		return(ans)
}


#given a triangular matrix, plot its heatmap, rotated 45 degrees (tip pointing downwards). if matrix is not triangular, values in the upper triangle are ignored. 
# x = R2.matrix. can be upper- or lower-triangular, so output of ld.script is fine.
# col_extremes specify colour range for the heatmap: c(col_of_low_LD, col_of_high_LD).
# 18.03.2015 --> v1.1

plot.heatmap <- function(x, col_extremes = c('khaki1', 'firebrick1')) {
  require('plotrix')
  #if matrix is not lower triangular: if it is upper triangular instead, transpose; else replace the upper triangle with NA
  if(all(!is.na(x[upper.tri(x)]))){
	if(all(is.na(x[lower.tri(x)]))) x <- t(x)
		else x[upper.tri(x)] <- NA }
  d <- nrow(x)
  p <- d/sqrt(2) #length of triangle sides needed to make hypotenuse of length d. then d/2 = is the height of the triangle 

  plot(NA, type = "n", xlim = c(0, d), ylim = c(0, d/2), asp = 1, axes = FALSE, xlab = "", ylab = "", bty = 'n')
  rasterImage(color.scale(x, extremes = col_extremes), xleft = d/2, xright = d/2 + p, ybottom = 0, ytop = p, interpolate = FALSE, angle = 45)
}


#given an LD matrix and number of SNPs per window produces an LD landscape plot (with indices on the x-axis) with the corresponding heatmap placed immediately underneath it
# 20.03.2015 --> v1.2
#add parameter d to determine length of the graphs? 

plot.landscape <- function(r2mat, n, col_extremes = c('khaki1', 'firebrick1'), type = 'l', lty = 2, col = 'dodgerblue', xlab = '', ylab = '', ...){
	ld <- ld.win(n = n, r2mat)
	layout(rbind(1, 2))
	par(oma = c(0, 0, 0, 0))
	par(mar = c(0, 2.5, 0, 0))
	plot(ld, type = type, lty = lty, col = col, xaxt = 'n', yaxt = 'n', bty = 'L', xlab = xlab, ylab = ylab, ...)
	axis(1, at = c(0, seq(ld)), labels = FALSE, tick = FALSE, lty = 1, col = 'grey30')
	axis(2, tick = TRUE, tck = -0.01, cex.axis = 0.85, las = 1)
	legend('topleft', legend = expression(R^2), lty = 2, col = col, bty = 'n')
	plot.heatmap(r2mat, col_extremes = col_extremes)
}
