setwd("******")
# Multi-core computing to save time
library(doParallel)
cl<-makeCluster(12)
registerDoParallel(cl)
source('wavelet-tool.R')

#
#### Figure-2 in result ####
# set the color palette
pal <- colorRampPalette(c('#003366', "#0072B2", "#56B4E9", '#FFFFFF', "#E69F00", "#D55E00", '#CC3300'), bias=1)(n=200)

tiff(filename = "Figure-2.tif", width = 5600, height = 2700, compression = "lzw", res = 300)
par(mfrow=c(2, 4))
for( i in c(60, 63, 66, 69) ){
  set.seed(123)
  raw.mat1 <- raw.mat2 <- matrix(NA, nrow = 101, ncol = 100)
  for(j in 1:100){
    raw.mat1[,j] <- sin(10*pi*seq(0,2,length.out = 101)+ circular::rvonmises(101, circular::circular(0), 10))
    raw.mat2[,j] <- sin(10*pi*seq(0,2,length.out = 101)+ circular::rvonmises(101, circular::circular(pi), 10))
  }
  data.mat <- cbind( raw.mat1[, 1:i] , raw.mat2[, (i+1):100] )
  w.time <- multi.wt(1:101, (data.mat))
  b.time <- colSums(abs(apply(w.time$z, c(1, 2), sum))) / colSums(apply(abs(w.time$z), c(1, 2), sum))
  b.time.boot <- b.boot(w.time, reps = 100, method = 'wmr')
  par(mar=c(4,3,3,3))
  plot(b.time, w.time$y, log='y', type = 'l', lwd=2,xaxs='i', xlim = c(-0.1,0.5), xlab = '', ylab = '', axes = F)
  axis(1, c(0,0.2,0.4) ) ; axis(4, c(2,5,10,20,40), line = -4.5 )
  lines(b.time.boot[,1], w.time$y, lty = 'dashed')
  lines(b.time.boot[,2], w.time$y, lty = 'dashed')
  title(xlab = expression(b[Time]), cex.lab = 1.5, font = 2, line = 2.5) ; mtext(4, text="Scale", line = -2)
  title(main = paste(i,'% and ', 100-i , '%'))
  mtext(paste('(',letters[(i-57)%/%3], ')',sep = ''), cex = 1.3, font = 2, line = 0.5, adj = -0.02)
}

i = 63 ; data.mat <- cbind(raw.mat1[,1:i] , raw.mat2[,(i+1):100] ) ; w.time <- multi.wt(1:101, (data.mat))
par(mar=c(4, 4, 2, 0.75))
image(1:101, 1:100, data.mat, col = pal, xlab = '', ylab = '', axes = F); box(col='black')
axis(1,  c(1, 21, 41, 61, 81, 101) ) ; axis(2, c(1, 21, 41, 61, 81, 100), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.3, font = 2, line = 2.8)
mtext('(e)', cex = 1.3, font = 2, line = 0.5, adj = -0.05)

Mb.time <- mb.com(w.time, 'wmr')
image(w.time$x, w.time$y, Mb.time, col = pal, log = 'y', xlab = '', ylab = '', axes = F); box(col = 'black')
axis(1,  c(1, 21, 41, 61, 81, 101) ) ; axis(2, c(2, 5, 10, 20, 40), las = 2)
title(xlab = 'Time', ylab = 'Scale', cex.lab = 1.3, font = 2, line = 2.8)
mtext('(f)', cex = 1.3, font = 2, line = 0.5, adj = -0.05)
time.boot <- mb.boot(w.time, reps = 100, method = 'wmr') ; time.line <- 1-abs(1-2*time.boot)
time.perm <- p.adjust(as.vector(time.line), method = "BY") ; dim(time.perm) = dim(time.line)
lines(min(w.time$x) + w.time$y, w.time$y, lty = 2, lwd = 1.5, col = "black")
lines(max(w.time$x) - w.time$y, w.time$y, lty = 2, lwd = 1.5, col = "black")
contour(w.time$x, w.time$y, time.line, levels = 0.01, lty = 3, add = TRUE, drawlabels = FALSE, axes=FALSE)
contour(w.time$x, w.time$y, time.perm, levels = 0.01, lwd = 2, add = TRUE, drawlabels = FALSE, axes=FALSE)

mat.t1 <- matrix(NA, 100, 101 )
for(i in 1:dim(w.time$z)[3]){  mat.t1[i,] <- w.time$z[, 51, i] }
image(1:101, 1:100, Mod(t(mat.t1)), col = pal, xlab = '', ylab = '', axes = F); box(col = 'black')
axis(1,  c(1, 21, 41, 61, 81, 101) ) ; axis(2, c(1, 21, 41, 61, 81, 100), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.3, font = 2, line = 2.8)
mtext('(g)', cex = 1.3, font = 2, line = 0.5, adj = -0.05)
image(1:101, 1:100, Re(t(mat.t1)), col = pal, xlab = '', ylab = '', axes = F); box(col='black')
axis(1,  c(1, 21, 41, 61, 81, 101) ) ; axis(2, c(1, 21, 41, 61, 81, 100), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.3, font = 2, line = 2.8)
mtext('(h)', cex = 1.3, font = 2, line = 0.5, adj = -0.05)
dev.off()

#
#### Figure-3 in result ####
sp.data <- openxlsx::read.xlsx('Fishing-data.xlsx', sheet=1)
sp.data.mat <- as.matrix(reshape2::dcast( sp.data, AREA.CODE ~ year , value.var = 'value')[, -1])
dat <- wsyn::cleandat(dat=as.matrix(sp.data.mat), times=c(1:60), clev=5)$cdat
sp.data.mat <- t(dat); colnames(sp.data.mat) <- unique(sp.data$AREA.CODE)

year <- unique(sp.data$year)
angle <- seq(min(sp.data$Latitude), max(sp.data$Latitude), length.out = 60)

sp.w.time <- multi.wt(1:nrow(sp.data.mat), (sp.data.mat))
sp.Mb.time <- mb.com(sp.w.time, 'wmr')
sp.time.boot <- mb.boot(sp.w.time, reps = 100, method = 'wmr') ; sp.time.line <- 1-abs(1-2*sp.time.boot)
sp.time.perm <- p.adjust(as.vector(sp.time.line), method = "BY") ; dim(sp.time.perm) = dim(sp.time.line)
sp.b.time <- colSums(abs(apply(sp.w.time$z, c(1, 2), sum))) / colSums(apply(abs(sp.w.time$z), c(1, 2), sum))
sp.b.time.boot <- b.boot(sp.w.time, reps = 100, method = 'wmr')

tiff(filename = "Figure-3.tif", width = 3600, height = 2000, compression = "lzw", res = 300)
par(mfrow = c(1,2), mar = c(4, 4, 0.75, 0.75), fig = c(0, 0.75, 0, 1))
image(year, sp.w.time$y, sp.Mb.time, log = 'y', col = pal, xlab = '', ylab = '', axes = F); box(col = 'black')
axis(1, c(1955, 1965, 1975, 1985, 1995, 2005, 2014) ) ; axis(2, c(2, 5, 13, 26), las = 2)
lines(min(year) + sp.w.time$y, sp.w.time$y, lty = 2, lwd = 1.5, col = "black")
lines(max(year) - sp.w.time$y, sp.w.time$y, lty = 2, lwd = 1.5, col = "black")
contour(year, sp.w.time$y, sp.time.line, levels = 0.01, lty = 3, add = TRUE, drawlabels = FALSE, axes=FALSE)
contour(year, sp.w.time$y, sp.time.perm, levels = 0.01, lwd = 2, add = TRUE, drawlabels = FALSE, axes=FALSE)
title(xlab = 'Time (year)', ylab = 'Scale (year)', cex.lab = 1.5, font = 2, line = 2.5)

par(mar = c(4, 0.75, 0.75, 0.75), fig = c(0.77, 1, 0, 1), new = T)
plot(sp.b.time, sp.w.time$y, log = 'y', type = 'l', lwd = 2, 
     yaxs = 'i', xaxs = 'i', xlim = c(0, 1), xlab = '', ylab = '', axes = F)
axis(1, c(0,0.2,0.4,0.6,0.8,1), labels = Vectorize(formatC)(c(0,0.2,0.4,0.6,0.8,1), digits=c(0,1,1,1,1,0)) )
lines(sp.b.time.boot[, 1], sp.w.time$y, lty = 'dashed')
lines(sp.b.time.boot[, 2], sp.w.time$y, lty = 'dashed')
title(xlab = expression(b[Time]), ylab = '', cex.lab = 1.5, font = 2, line = 3)
dev.off()

#
#### Figure-4 in result ####
# Extend raw data and add another scale
obj <- list(x = seq(nrow(sp.data.mat)), y = seq(ncol(sp.data.mat)), z = (sp.data.mat))
grid.list<- list( x = seq(1, nrow(sp.data.mat), , nrow(sp.data.mat)), y= seq(1, ncol(sp.data.mat), , nrow(sp.data.mat)))
newmat <- fields::interp.surface.grid(obj, grid.list )

scales <- 1:nrow(sp.data.mat)
set.seed(12)
raw.mat1 <- raw.mat11 <- raw.mat2 <- matrix(NA, nrow = length(scales), ncol=length(scales))
for(i in 1:(length(scales))){
  raw.mat1[i,] <- sin(10*pi*seq(0, 2, length.out = length(scales))+ circular::rvonmises(length(scales), circular::circular(0), 15))
  raw.mat11[i[i%%2==0], ] <- sin(4*pi*seq(0, 2, length.out = length(scales))+ circular::rvonmises(length(scales), circular::circular(0), 15))
  raw.mat11[i[i%%2==1], ] <- sin(4*pi*seq(0, 2, length.out = length(scales))+ circular::rvonmises(length(scales), circular::circular(pi), 15))
  raw.mat2[, i] <- cos(circular::rvonmises(length(scales), circular::circular(sample(0:360, 1)*pi/180), 2))
}
data.mat <- newmat$z + rbind(raw.mat1[1:30, ] , raw.mat11[1:30, ])*2 + raw.mat2

w.time <- multi.wt(scales, (data.mat)) ; w.space <- multi.wt(scales, t(data.mat))
Mb.time <- mb.com(w.time, 'wmr') ; Mb.space <- mb.com(w.space, 'wmr')
b.time <- colSums(abs(apply(w.time$z, c(1, 2), sum))) / colSums(apply(abs(w.time$z), c(1, 2), sum))
b.space <- colSums(abs(apply(w.space$z, c(1, 2), sum))) / colSums(apply(abs(w.space$z), c(1, 2), sum))

# Figure-4-1:
tiff(filename = "Figure-4-1.tif", width = 3000, height = 3300, compression = "lzw", res = 300)
par(mar=c(4, 4, 1.5, 0.75), fig=c(0, 0.75, 0.5, 1))
image(year, w.time$y, Mb.time, log = 'y', col = pal, xlab = '', ylab = '', axes = F); box(col = 'black')
axis(1, c(1955, 1965, 1975, 1985, 1995, 2005, 2014) ) ; axis(2, c(2, 5, 13, 26), las = 2)
lines(min(year) + w.time$y, w.time$y, lty = 2, lwd = 1.5, col = "black")
lines(max(year) - w.time$y, w.time$y, lty = 2, lwd = 1.5, col = "black")
title(xlab = 'Time (year)', ylab = 'Scale (year)', cex.lab = 1.5, font = 2, line = 2.5)
mtext('(a)',cex = 1.3, font = 2, line = 0.5, adj = -0.1)

par(mar=c(4, 0.5, 1.5, 0.75), fig=c(0.77, 1, 0.5, 1), new=T)
plot(b.time, w.time$y, log = 'y', type = 'l', lwd = 2, 
     yaxs = 'i',xaxs = 'i', xlim = c(0,1), xlab = '', ylab = '', axes = F)
axis(1, c(0,0.2,0.4,0.6,0.8,1), labels = Vectorize(formatC)(c(0,0.2,0.4,0.6,0.8,1), digits=c(0,1,1,1,1,0)) )
title(xlab = expression(b[Time]), ylab = '', cex.lab = 1.5, font = 2, line = 3)

par(mar=c(4,4,1.5,0.75), fig=c(0,0.75,0,0.5), new=T)
image(angle, w.space$y, Mb.space, log='y', col=pal, xlab='', ylab='', axes = F); box(col='black')
axis(1, c(-30,-15,0,15,30,45), labels = c('30S','15S','0','15N','30N','45N') )
axis(2, c(3, 6, 16), las = 2)
lines(min(angle) + w.space$y, w.space$y, lty = 2, lwd = 1.5, col = "black")
lines(max(angle) - w.space$y, w.space$y, lty = 2, lwd = 1.5, col = "black")
title(xlab = 'Space (°)', ylab = 'Scale (°)', cex.lab = 1.5, font = 2, line = 2.5) 
mtext('(b)',cex = 1.3, font = 2, line = 0.5, adj = -0.1)

par(mar = c(4, 0.5, 1.5, 0.75), fig = c(0.77, 1, 0, 0.5), new = T)
plot(b.space, w.space$y, log = 'y', type = 'l', lwd = 2, 
     yaxs = 'i',xaxs = 'i', xlim = c(0,1), xlab = '', ylab = '', axes = F)
axis(1, c(0,0.2,0.4,0.6,0.8,1), labels = Vectorize(formatC)(c(0,0.2,0.4,0.6,0.8,1), digits=c(0,1,1,1,1,0)) )
title(xlab = expression(b[Space]), ylab = '', cex.lab = 1.5, font = 2, line = 3)
dev.off()

# Figure-4-2:
mat <- outer(maxmin(b.space), maxmin(b.time), '*')
tiff(filename = "Figure-4-2.tif", width = 2100, height = 1600, compression = "lzw", res = 300)
par(mar = c(4, 4.5, 1.5, 0.75))
image(w.time$y, w.space$y, t(mat), xlab = '', ylab = '', col = pal, log = 'xy', axes = F); box(col = 'black')
axis(1, c(2, 5, 13, 26) ) ; axis(2, c(3, 6, 16), las = 2)
rect(10.5, 4, 16, 8.5, border = "green", lwd = 3) ; rect(10.5, 11, 16, 23, border = "red", lwd = 3)
title(xlab = expression(b[Time]), ylab = expression(b[Space]), cex.lab = 1.5, font = 2, line = 2.5)
mtext('(c)',cex = 1.3, font = 2, line = 0.5, adj = -0.1)
dev.off()

#
#### Figure-5 in result ####
mat.t1 <- mat.s1 <- mat.s2 <- matrix(NA, length(scales), length(scales) )
for(i in 1:dim(w.time$z)[3]){  mat.t1[i, ] <- w.time$z[,42,i] }
for(i in 1:dim(w.space$z)[3]){  mat.s1[,i ] <- w.space$z[,25,i] ; mat.s2[,i ] <- w.space$z[,48,i] }
tiff(filename = "Figure-5.tif", width = 3600, height = 1200, compression = "lzw", res = 300)
par(mfrow = c(1,3), mar = c(4, 4, 2, 0.75))
image(year, angle, (data.mat), col = pal, axes = FALSE, xlab = '', ylab = ''); box(col = 'black')
title(xlab = 'Time (year)', ylab = 'Space (°)', cex.lab = 1.5, cex.main = 1.5, font = 2, line = 2.8)
axis(1,  c(1955, 1965, 1975, 1985, 1995, 2005, 2014) )
axis(2, c(-30, -15, 0, 15, 30, 45), labels = c('30S','15S','0','15N','30N','45N'), las = 2)
mtext('(a)',cex = 1.3, font = 2, line = 0.5, adj = -0.1)
image(year, angle, sqrt(Mod(t(mat.t1))*Mod(t(mat.s1))), col = pal, axes = FALSE, xlab='', ylab=''); box(col='black')
title(xlab = 'Time (year)', ylab = 'Space (°)', cex.lab = 1.5, cex.main = 1.5, font = 2, line = 2.8)
axis(1,  c(1955, 1965, 1975, 1985, 1995, 2005, 2014) )
axis(2, c(-30, -15, 0, 15, 30, 45), labels = c('30S','15S','0','15N','30N','45N'), las = 2)
mtext('(b)',cex = 1.3, font = 2, line = 0.5, adj = -0.1)
image(year, angle, sqrt(Mod(t(mat.t1))*Mod(t(mat.s2))), col=pal, axes = FALSE, xlab='', ylab=''); box(col='black')
title(xlab = 'Time (year)', ylab = 'Space (°)', cex.lab = 1.5, cex.main = 1.5, font = 2, line = 2.8)
axis(1,  c(1955, 1965, 1975, 1985, 1995, 2005, 2014) )
axis(2, c(-30, -15, 0, 15, 30, 45), labels = c('30S','15S','0','15N','30N','45N'), las = 2)
mtext('(c)',cex = 1.3, font = 2, line = 0.5, adj = -0.1)
dev.off()

#
#### Figure-6 a-d in result ####
times <- 1:600; spaces <- 1:400
set.seed(1234)
raw.mat1 <- raw.mat11 <- matrix(NA, nrow = length(spaces), ncol=length(times))
raw.mat2 <- raw.mat22 <- matrix(NA, nrow = length(times), ncol=length(times))
for(j in 1:(length(times))){
  raw.mat1[,j] <- sin(10*pi*seq(0,2,length.out = length(spaces))+ circular::rvonmises(length(spaces), circular::circular(0), 5))
  raw.mat11[,j] <- sin(40*pi*seq(0,2,length.out = length(spaces))+ circular::rvonmises(length(spaces), circular::circular(0), 5))
  raw.mat2[j,] <- sin(10*pi*seq(0,2,length.out = length(times))+ circular::rvonmises(length(times), circular::circular(0), 5))
  raw.mat22[j,] <- sin(40*pi*seq(0,2,length.out = length(times))+ circular::rvonmises(length(times), circular::circular(0), 5))
}
data.mat <- rbind( cbind((raw.mat11[1:200,1:300] + raw.mat22[1:200,1:300]) , (raw.mat11[1:200,301:600] + raw.mat2[1:200,301:600]) ),
                   cbind((raw.mat1[201:400,1:300] + raw.mat22[201:400,1:300]) , (raw.mat1[201:400,301:600] + raw.mat2[201:400,301:600]) ) )

w.time <- multi.wt(times, t(data.mat)) ; w.space <- multi.wt(spaces, data.mat)
b.time <- colSums(abs(apply(w.time$z, c(1, 2), sum)))/colSums(apply(abs(w.time$z), c(1, 2), sum))
b.space <- colSums(abs(apply(w.space$z, c(1, 2), sum)))/colSums(apply(abs(w.space$z), c(1, 2), sum))
mat <- outer((b.space), (b.time), '*')
#
tiff(filename = "Figure-6-1.tif", width = 3200, height = 4000, compression = "lzw", res = 300)
par(mar=c(4, 4.3, 2, 0.75), fig=c(0, 1, 0.4, 1))
image(1:600, 1:400, t(data.mat), col = pal, xlab = '', ylab = '', axes = F); box(col = 'black')
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.5, font = 2, line = 2.5)
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
mtext('(a)',cex = 1.5, font = 2, line = 0.5, adj = -0.08)

par(mar=c(4, 1, 1.5, 0.75), fig=c(0, 0.2, 0, 0.4), new = T)
plot(b.time, w.time$y, log = 'y', type = 'l', lwd = 2, 
     yaxs = 'i',xaxs = 'i', xlim = c(0,1), xlab = '', ylab = '', axes = F)
axis(1, c(0, 0.5, 1), labels = Vectorize(formatC)(c(0.0,0.5,1.0), digits=c(0,1,0)) )
title(xlab = expression(b[Time]), ylab = '', cex.lab = 1.5, font = 2, line = 3)
mtext('(b)',cex = 1.5, font = 2, line = 0.5, adj = -0.05)

par(mar=c(4, 1, 1.5, 0.75), fig=c(0.2, 0.4, 0, 0.4), new = T)
plot(b.space, w.space$y, log = 'y', type = 'l', lwd = 2,
     yaxs = 'i',xaxs = 'i', xlim = c(0,1), xlab = '', ylab = '', axes = F)
axis(1, c(0, 0.5, 1), labels = Vectorize(formatC)(c(0.0,0.5,1.0), digits=c(0,1,0)) )
title(xlab = expression(b[Space]), ylab = '', cex.lab = 1.5, font = 2, line = 3)

par(mar=c(4, 4.5, 1.5, 0.75), fig=c(0.4, 1, 0, 0.4), new = T)
image(w.time$y, w.space$y, t(mat), xlab = '', ylab = '', col = pal, log = 'xy', axes = F); box(col = 'black')
axis(1, c(2, 5, 15, 60, 240) ) ; axis(2, c(2, 10, 40, 160), las = 2)
title(xlab = expression(b[Time]), ylab = expression(b[Space]), cex.lab = 1.5, font = 2, line = 2.5)
mtext('(c)',cex = 1.5, font = 2, line = 0.5, adj = -0.1)
dev.off()

# 
mat.t1 <- mat.t2 <- mat.s1 <- mat.s2 <- matrix(NA, length(spaces), length(times) )
for(i in 1:dim(w.time$z)[3]){  mat.t1[i, ] <- w.time$z[, 241, i] ; mat.t2[i, ] <- w.time$z[, 408, i] }
for(i in 1:dim(w.space$z)[3]){  mat.s1[, i] <- w.space$z[, 140, i] ; mat.s2[, i] <- w.space$z[, 261, i] }
tiff(filename = "Figure-6-2.tif", width = 4000, height = 3000, compression = "lzw", res = 300)
par(mfcol = c(2, 2), mar = c(4, 4, 3, 0.75))
image(times, spaces, sqrt(Mod(t(mat.t1))*Mod(t(mat.s2))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
mtext('(d)',cex = 1.5, font = 2, line = 1, adj = -0.1)
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Short time scale and large space scale', cex.main = 1.5, font = 2, line = 0.6)
image(times, spaces, sqrt(Mod(t(mat.t1))*Mod(t(mat.s1))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Short time scale and small space scale', cex.main = 1.5, font = 2, line = 0.6)
image(times, spaces, sqrt(Mod(t(mat.t2))*Mod(t(mat.s2))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Long time scale and large space scale', cex.main =1.5, font = 2, line = 0.6)
image(times, spaces, sqrt(Mod(t(mat.t2))*Mod(t(mat.s1))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Long time scale and small space scale', cex.main = 1.5, font = 2, line = 0.6)
dev.off()

#
#### Figure-6 e-h in result ####
times <- 1:600; spaces <- 1:400
set.seed(1234)
raw.mat1 <- raw.mat11 <- matrix(NA, nrow = length(spaces), ncol=length(times))
raw.mat2 <- raw.mat22 <- raw.mat3 <- matrix(NA, nrow = length(times), ncol = length(times))
for(j in 1:(length(times))){
  raw.mat1[, j] <- sin(10*pi*seq(0, 2, length.out = length(spaces))+ circular::rvonmises(length(spaces), circular::circular(0), 20))
  raw.mat11[, j[j%%2==0]] <- sin(40*pi*seq(0, 2, length.out = length(spaces))+ circular::rvonmises(length(spaces), circular::circular(0), 20))
  raw.mat11[, j[j%%2==1]] <- sin(40*pi*seq(0, 2, length.out = length(spaces))+ circular::rvonmises(length(spaces), circular::circular(pi), 20))
  raw.mat2[j, ] <- sin(10*pi*seq(0, 2, length.out = length(times))+ circular::rvonmises(length(times), circular::circular(0), 20))
  raw.mat22[j[j%%2==0], ] <- sin(40*pi*seq(0, 2, length.out = length(times))+ circular::rvonmises(length(times), circular::circular(0), 20))
  raw.mat22[j[j%%2==1], ] <- sin(40*pi*seq(0, 2, length.out = length(times))+ circular::rvonmises(length(times), circular::circular(pi), 20))
  raw.mat3[j, ] <- cos(circular::rvonmises(length(times), circular::circular(sample(0:360, 1)*pi/180), 5))
}
data.mat <- rbind( cbind((6*raw.mat11[1:200, 1:300] + 6*raw.mat22[1:200, 1:300]) , (6*raw.mat11[1:200, 301:600] + raw.mat2[1:200, 301:600]) ),
                   cbind((raw.mat1[201:400, 1:300] + 6*raw.mat22[201:400, 1:300]) , (raw.mat1[201:400, 301:600] + raw.mat2[201:400, 301:600]) ) ) +
  raw.mat3[1:length(spaces), 1:length(times)] + t(raw.mat3[1:length(times), 1:length(spaces)])

w.time <- multi.wt(times, t(data.mat)) ; w.space <- multi.wt(spaces, data.mat)
b.time <- colSums(abs(apply(w.time$z, c(1, 2), sum))) / colSums(apply(abs(w.time$z), c(1, 2), sum))
b.space <- colSums(abs(apply(w.space$z, c(1, 2), sum))) / colSums(apply(abs(w.space$z), c(1, 2), sum))
mat <- outer(maxmin(b.space), maxmin(b.time), '*')
# 
tiff(filename = "Figure-7-1.tif", width = 3200, height = 4000, compression = "lzw", res = 300)
par(mar=c(4, 4.3, 2, 0.75), fig=c(0, 1, 0.4, 1))
image(1:600, 1:400, t(data.mat), col = pal, xlab = '', ylab = '', axes = F); box(col = 'black')
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.5, font = 2, line = 2.5)
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
mtext('(a)',cex = 1.5, font = 2, line = 0.5, adj = -0.08)

par(mar=c(4, 1, 1.5, 0.75), fig=c(0, 0.2, 0, 0.4), new = T)
plot(b.time, w.time$y, log = 'y', type = 'l', lwd = 2, 
     yaxs = 'i',xaxs = 'i', xlim = c(0,1), xlab = '', ylab = '', axes = F)
axis(1, c(0, 0.5, 1), labels = Vectorize(formatC)(c(0.0,0.5,1.0), digits=c(0,1,0)) )
title(xlab = expression(b[Time]), ylab = '', cex.lab = 1.5, font = 2, line = 3)
mtext('(b)',cex = 1.5, font = 2, line = 0.5, adj = -0.05)

par(mar=c(4, 1, 1.5, 0.75), fig=c(0.2, 0.4, 0, 0.4), new = T)
plot(b.space, w.space$y, log = 'y', type = 'l', lwd = 2,
     yaxs = 'i',xaxs = 'i', xlim = c(0,1), xlab = '', ylab = '', axes = F)
axis(1, c(0, 0.5, 1), labels = Vectorize(formatC)(c(0.0,0.5,1.0), digits=c(0,1,0)) )
title(xlab = expression(b[Space]), ylab = '', cex.lab = 1.5, font = 2, line = 3)

par(mar=c(4, 4.5, 1.5, 0.75), fig=c(0.4, 1, 0, 0.4), new = T)
image(w.time$y, w.space$y, t(mat), xlab = '', ylab = '', col = pal, log = 'xy', axes = F); box(col = 'black')
axis(1, c(2, 5, 15, 60, 240) ) ; axis(2, c(2, 10, 40, 160), las = 2)
title(xlab = expression(b[Time]), ylab = expression(b[Space]), cex.lab = 1.5, font = 2, line = 2.5)
mtext('(c)',cex = 1.5, font = 2, line = 0.5, adj = -0.1)
dev.off()

#
mat.t1 <- mat.t2 <- mat.s1 <- mat.s2 <- matrix(NA, length(spaces), length(times) )
for(i in 1:dim(w.time$z)[3]){  mat.t1[i, ] <- w.time$z[, 241, i] ; mat.t2[i, ] <- w.time$z[, 408, i] }
for(i in 1:dim(w.space$z)[3]){  mat.s1[, i] <- w.space$z[, 140, i] ; mat.s2[, i] <- w.space$z[, 261, i] }
tiff(filename = "Figure-7-2.tif", width = 4000, height = 3000, compression = "lzw", res = 300)
par(mfcol = c(2, 2), mar = c(4, 4, 3, 0.75))
image(times, spaces, sqrt(Mod(t(mat.t1))*Mod(t(mat.s2))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
mtext('(d)',cex = 1.5, font = 2, line = 1, adj = -0.1)
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Short time scale and large space scale', cex.main = 1.5, font = 2, line = 0.6)
image(times, spaces, sqrt(Mod(t(mat.t1))*Mod(t(mat.s1))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Short time scale and small space scale', cex.main = 1.5, font = 2, line = 0.6)
image(times, spaces, sqrt(Mod(t(mat.t2))*Mod(t(mat.s2))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Long time scale and large space scale', cex.main =1.5, font = 2, line = 0.6)
image(times, spaces, sqrt(Mod(t(mat.t2))*Mod(t(mat.s1))), 
      xlab = '', ylab = '', col = pal, axes = F); box(col = 'black')
axis(1,  c(1, 100, 200, 300, 400, 500, 600) ) ; axis(2, c(1, 100, 200, 300, 400), las = 2)
title(xlab = 'Time', ylab = 'Space', cex.lab = 1.2, font = 2, line = 2.5)
title(main = 'Long time scale and small space scale', cex.main = 1.5, font = 2, line = 0.6)
dev.off()

#
#### End ####
#
stopCluster(cl)


