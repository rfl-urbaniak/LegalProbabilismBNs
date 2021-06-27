
Rafal <- function(px,pyx,pynx){
  pnx <- 1-px
  py <- pyx*px + pynx * pnx 
  pny <- 1-py
  pnyx <- 1-pyx
  pnynx <- 1-pynx
  Z(pyx,py) * px+ Z(pynx,py) * pnx + Z(pnyx,pny) * px + Z(pnynx,pny) * pnx
}

Alicja <-  function(px,pyx,pynx){
  pnx <- 1-px
  py <- pyx*px + pynx *  pnx
  pny <- 1-py
  pnyx <- 1-pyx
  pnynx <- 1-pynx
  Z(pyx,py) * px * py + Z(pynx,py) * pnx * py + Z(pnyx,pny) * px * pny + Z(pnynx,pny) * pnx * pny
}

Joint <-  function(px,pyx,pynx){
  pnx <- 1-px
  py <- pyx*px + pynx *  pnx
  pny <- 1-py
  pnyx <- 1-pyx
  pnynx <- 1-pynx
  XY <- pyx * px
  XnY <- pnyx * px
  nXY <-pynx * pnx
  nXnY <- pnynx * pnx
  Z(pyx,py) * XY  + Z(pynx,py) * nXY + Z(pnyx,pny) * XnY + Z(pnynx,pny) * nXnY
}

pyx <- seq(0,1,by=0.01)
pynx <- seq(0,1,by=0.01)
options <- expand.grid(pyx,pynx)
colnames(options) <- c("pyx","pynx")

options$r1 <- Rafal(.1,options$pyx,options$pynx)
options$a1 <- Alicja(.1,options$pyx,options$pynx)
options$j1 <- Joint(.1,options$pyx,options$pynx)




rplot1 <- scatter3D(options$pyx,options$pynx,options$r1,pch=19,cex=.2,byt="g",alpha=0.6,theta=30, phi=20,ticktype="simple",xlab="pyx", ylab="pynx",zlab="rafal",main="Rafal,  prior x =.1",colvar=NULL)


aplot1 <- scatter3D(options$pyx,options$pynx,options$a1,pch=19,cex=.2,byt="g",alpha=.6,theta=40, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="alicja",  main="Alicja, prior x =.1",colvar=NULL)

jplot1 <- scatter3D(options$pyx,options$pynx,options$j1,pch=19,cex=.2,byt="g",alpha=.6,theta=40, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="joint",
                    main="Joint, prior x =.1",colvar=NULL)


options$r2 <- Rafal(.5,options$pyx,options$pynx)
options$a2 <- Alicja(.5,options$pyx,options$pynx)
options$j2 <- Joint(.5,options$pyx,options$pynx)



rplot2 <- scatter3D(options$pyx,options$pynx,options$r2,pch=19,cex=.2,byt="g",alpha=.6,theta=40, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="rafal",
                    main="Rafal, prior x =.5",colvar=NULL)

aplot2 <- scatter3D(options$pyx,options$pynx,options$a2,pch=19,cex=.2,byt="g",alpha=.6,theta=20, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="alicja",
                    main="Alicja, prior x =.5",colvar=NULL)

jplot2 <- scatter3D(options$pyx,options$pynx,options$j2,pch=19,cex=.2,byt="g",alpha=.6,theta=25, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="joint",
                    main="Joint, prior x =.5",colvar=NULL)


#aplot2


options$r3 <- Rafal(.9,options$pyx,options$pynx)
options$a3 <- Alicja(.9,options$pyx,options$pynx)
options$j3 <- Joint(.9,options$pyx,options$pynx)


rplot3 <- scatter3D(options$pyx,options$pynx,options$r2,pch=19,cex=.2,byt="g",alpha=.6,theta=40, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="rafal",
                    main="Rafal, prior x =.9",colvar=NULL)

aplot3 <- scatter3D(options$pyx,options$pynx,options$a2,pch=19,cex=.2,byt="g",alpha=.6,theta=20, phi=30,ticktype="simple",xlab="pyx",ylab="pynx",zlab="alicja",
                    main="Alicja, prior x =.9",colvar=NULL)

jplot3 <- scatter3D(options$pyx,options$pynx,options$j2,pch=19,cex=.2,byt="g",alpha=.6,theta=15, phi=40,ticktype="simple",xlab="pyx",ylab="pynx",zlab="joint",
                    main="Joint, prior x =.9",colvar=NULL)


#jplot1
#jplot2 
#jplot3

