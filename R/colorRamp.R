
colorRampPalette<-function(palette, bias=1,method=c("linear","spline")){
	
	if (bias<=0) stop("bias must be positive")
	coord<-as.data.frame(t(col2rgb(palette))/255)
	x<-seq(0,1,length=length(palette))
	if(match.arg(method)=="spline"){
	  r<-splinefun(x,coord$red)
	  g<-splinefun(x,coord$green)
	  b<-splinefun(x,coord$blue)
	}else{
	  r<-approxfun(x,coord$red)
	  g<-approxfun(x,coord$green)
	  b<-approxfun(x,coord$blue)
	}

       function(n){
	   x<-seq(0,1,length=n)^bias
	   rgb(r(x),g(x),b(x))
	}

}

colorRamp<-function(palette, bias=1,method=c("linear","spline")){
	
	coord<-as.data.frame(t(col2rgb(palette))/255)
	x<-seq(0,1,length=length(palette))
	if(match.arg(method)=="spline"){
	  r<-splinefun(x,coord$red)
	  g<-splinefun(x,coord$green)
	  b<-splinefun(x,coord$blue)
	}else{
	  r<-approxfun(x,coord$red)
	  g<-approxfun(x,coord$green)
	  b<-approxfun(x,coord$blue)
	}

        function(x,min=0,max=1){
	   if(any(x<min | x>max)) 
		stop("out of range")
	   z<-((x-min)/(max-min))^bias
	   rgb(r(z),g(z),b(z))
        }

}
