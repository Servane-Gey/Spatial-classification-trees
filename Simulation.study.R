# Simulation experiments from the article 'Spatial classification trees'
library(spatcart)
library(tree)
library(spatstat)
library(ggplot2)

set.seed(12)

# Data simulation
n = 1000
### Chess data set
nom="Chess"
chess = damier(n, h=0.45, model="Poisson")

#### Locally repulsive data set
# nom = "Repulsion"
# r0 = 0.05
# repuls = repulsion(n,r0)

donnees = chess$data
# donnees = repuls

# Define data as a point process with spatstat
ypp = ppp(x=donnees$x1, y=donnees$x2, marks = donnees$label,
          range(donnees$x1), range(donnees$x2))

major =  names(which.max(intensity(ypp)))
minor = names(which.min(intensity(ypp)))

K = Kcross(ypp, major, minor, correction="none")

if(nom != "Repulsion"){
  r = K$r[which.min(K$un-K$theo)]
} else {
  r = r0-r0/10
}

K01 = data.frame(cbind(K$r, K$un - K$theo))
colnames(K01) = c("r","difference")

DiffK = ggplot(K01, aes(x=r, y=difference))+geom_line()+
  geom_vline(xintercept = r, colour = "blue", linetype = "dashed")+
  xlab("r")+
  ylab("Kest - Ktheo")+
  ggtitle("Difference between estimated and theoretical K")
DiffK

###################################################
#              Trees' construction                #
###################################################

# Set stopping rule parameters in maximal tree
minsplit = 50
minleaf = 25

##### SpatCART
# Pruning with Gini deviance
tss_dev = spatcart(ypp, r, minsplit=minsplit,minleaf=minleaf, graph = FALSE)

# Pruning with misclassification rate
tss_mcr = spatcart(ypp, r, method="misclass", minsplit=minsplit,minleaf=minleaf, graph = FALSE)

##### CART
donnees2=data.frame(ypp)
t = tree(marks~.,data=donnees2,split="gini",model=T, minsize=minsplit,mincut=minleaf)

# Pruning with Gini deviance
seq_dev = prune.tree(t)

nbleav = choice.tree(seq_dev)
t_dev.max = prune.tree(t, best=nbleav$maxgap)
t_dev.min = prune.tree(t, best = nbleav$plateau)

# Pruning with misclassification rate
seq_mcr = prune.misclass(t)

nbleav = choice.tree(seq_mcr)
t_mcr.max = prune.misclass(t, best=nbleav$maxgap)
t_mcr.min = prune.misclass(t, best=nbleav$plateau)

###################################################
#           Graphical representations             #
###################################################

donnees2=data.frame(ypp)

###########################
# Class probability trees
#

tss = tss_dev
t1 = t_dev.max
t2 = t_dev.min
seq = seq_dev

#############################
# Dimension versus complexity

#### SpatCART
nbleav.maxgap = sum(tss$opt.tree.max$frame$var=="<leaf>")
nbleav.plateau = sum(tss$opt.tree.min$frame$var=="<leaf>")
alpha.maxgap = tss$pruned.seq$k[tss$pruned.seq$size==nbleav.maxgap]
alpha.plateau = tss$pruned.seq$k[tss$pruned.seq$size==nbleav.plateau]

plot(stepfun(tss$pruned.seq$k[-1],tss$pruned.seq$size), main=c(paste(nom), "SpatCART Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(-30,alpha.maxgap), c(nbleav.maxgap, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(alpha.plateau,alpha.plateau), c(0, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")
lines(c(-30,alpha.plateau), c(nbleav.plateau, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")

#### CART
nbleav.maxgap = sum(t1$frame$var=="<leaf>")
nbleav.plateau = sum(t2$frame$var=="<leaf>")
alpha.maxgap = seq$k[seq$size==nbleav.maxgap]
alpha.plateau = seq$k[seq$size==nbleav.plateau]

plot(stepfun(seq$k[-1],seq$size), main=c(paste(nom), "CART Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(-30,alpha.maxgap), c(nbleav.maxgap, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(alpha.plateau,alpha.plateau), c(0, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")
lines(c(-30,alpha.plateau), c(nbleav.plateau, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")

##################
# Partitions

#### SpatCART
# Optimal SpatCART trees
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART largest optimal tree"), cex.main = 1.5)
partition.tree.corrected(tss$opt.tree.max, add=T, ordvar=c("x","y"))

plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART smallest optimal tree"), cex.main = 1.5)
partition.tree.corrected(tss$opt.tree.min, add=T, ordvar=c("x","y"))

# If trees are the same, change the figure's title
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART optimal tree"), cex.main = 1.5)
partition.tree.corrected(tss$opt.tree.min, add=T, ordvar=c("x","y"))

# SpatCART maximal tree
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART maximal tree"), cex.main = 1.5)
partition.tree.corrected(tss$max.tree, add=T, ordvar=c("x","y"))

#### CART
# Optimal CART trees
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", cex.lab = 2, main=c(paste(nom), " CART largest optimal tree"), cex.main = 1.5)
partition.tree.corrected(t1, add=T, ordvar = c("x","y"))

plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", cex.lab = 2, main=c(paste(nom), " CART smallest optimal tree"), cex.main = 1.5)
partition.tree.corrected(t2, add=T, ordvar = c("x","y"))

# If trees are the same, change the figure's title
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", cex.lab = 2, main=c(paste(nom), " CART optimal tree"), cex.main = 1.5)
partition.tree.corrected(t2, add=T, ordvar = c("x","y"))

# Maximal CART tree
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", cex.lab = 2, main=c(paste(nom), " CART maximal tree"), cex.main = 1.5)
partition.tree.corrected(t, add=T, ordvar = c("x","y"))

######################
# Classification trees
#

tss = tss_mcr
t1 = t_mcr.max
t2 = t_mcr.min
seq = seq_mcr

#############################
# Dimension versus complexity

#### SpatCART
nbleav.maxgap = sum(tss$opt.tree.max$frame$var=="<leaf>")
nbleav.plateau = sum(tss$opt.tree.min$frame$var=="<leaf>")
alpha.maxgap = tss$pruned.seq$k[tss$pruned.seq$size==nbleav.maxgap]
alpha.plateau = tss$pruned.seq$k[tss$pruned.seq$size==nbleav.plateau]

plot(stepfun(tss$pruned.seq$k[-1],tss$pruned.seq$size), main=c(paste(nom), "SpatCART Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(-40,alpha.maxgap), c(nbleav.maxgap, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(alpha.plateau,alpha.plateau), c(0, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")
lines(c(-40,alpha.plateau), c(nbleav.plateau, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")

#### CART
nbleav.maxgap = sum(t1$frame$var=="<leaf>")
nbleav.plateau = sum(t2$frame$var=="<leaf>")
alpha.maxgap = seq$k[seq$size==nbleav.maxgap]
alpha.plateau = seq$k[seq$size==nbleav.plateau]

plot(stepfun(seq$k[-1],seq$size), main=c(paste(nom), "CART Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(-30,alpha.maxgap), c(nbleav.maxgap, nbleav.maxgap), type = "b", pch = 2, lty="dashed", col="blue")
lines(c(alpha.plateau,alpha.plateau), c(0, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")
lines(c(-30,alpha.plateau), c(nbleav.plateau, nbleav.plateau), type = "b", pch = 8, lty="dotdash", col="red")

#############
# Partitions

### SpatCART
# SpatCART optimal tree
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART max optimal tree"), cex.main = 1.5)
partition.tree.corrected(tss$opt.tree.max, add=T, ordvar=c("x","y"))

plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART min optimal tree"), cex.main = 1.5)
partition.tree.corrected(tss$opt.tree.min, add=T, ordvar=c("x","y"))

# If trees are the same, change the figure's title
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART optimal tree"), cex.main = 1.5)
partition.tree.corrected(tss$opt.tree.min, add=T, ordvar=c("x","y"))

# SpatCART maximal tree
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " SpatCART maximal tree"), cex.main = 1.5)
partition.tree.corrected(tss$max.tree, add=T, ordvar=c("x","y"))

### CART
# Optimal CART tree
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " CART optimal classification tree"), cex.main = 1.5)
partition.tree.corrected(t1, add=T, ordvar = c("x","y"))

plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " CART min optimal tree"), cex.main = 1.5)
partition.tree.corrected(t2, add=T, ordvar = c("x","y"))

# If trees are the same, change the figure's title
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", cex.lab = 2, main=c(paste(nom), " CART optimal tree"), cex.main = 1.5)
partition.tree.corrected(t2, add=T, ordvar = c("x","y"))

# Maximal CART tree
plot(ypp$x,ypp$y, pch=19, cex=0.4, col=c("blue","red")[ypp$marks], xlab="", ylab="", main=c(paste(nom), " CART maximal tree"), cex.main = 1.5)
partition.tree.corrected(t, add=T, ordvar = c("x","y"))


