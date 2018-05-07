# PCA PLOTS 3D
## in your local machine
```{r}
library(rgl)

====================
  TAM PCAplots

## PCAplot TAM_PostSepsis_TAM.naive
dat.pca <- prcomp(t(xyx))
dat.loads <- dat.pca$x[,1:3]
plot3d(dat.loads, size = 15, col = c("green", "green", "green", "green", "orange", "orange", "orange"))
rgl.postscript("plot3d.tam.pdf", "pdf")

## PCAplotTAM_PostSepsis_M1
dat.pca.m1 <- prcomp(t(xyx.m1))
dat.loads.m1 <- dat.pca.m1$x[,1:3]
plot3d(dat.loads.m1, size = 15, col = c("green", "green", "green", "green", "black", "black", "black", "black"))
rgl.postscript("plot3d.m1.pdf", "pdf")

## PCAplot TAM_PostSepsis_M2
dat.pca.m2 <- prcomp(t(xyx.m2))
dat.loads.m2 <- dat.pca.m2$x[,1:3]
plot3d(dat.loads.m2, size = 15, col = c("green", "green", "green", "green", "red", "red", "red", "red"))
rgl.postscript("plot3d.m2.pdf", "pdf")

====
  Movie from PCAplot
## Only for on my own computer - rgl code
plot3d(dat.loads, col = c("green", "green", "green", "green", "red", "red", "red", "red"), type='s')
rgl.clear(type="lights")
rgl.light(-45, 20, ambient="black", diffuse="#dddddd", specular="white")
rgl.light(60, 30, ambient="#dddddd", diffuse="#dddddd", specular="black")
## before you animate your movie, make the window big
movie3d(spin3d(), fps=60, duration=10, dir="/home/thais/",type = "gif")
## Then not in R, look in /tmp/ for a directory that starts with R
## enter that directory and do this:
ffmpeg -f image2 -r 30 -pattern_type glob -i '*.png' movie.mp4
cp movie.mp4 ~/
