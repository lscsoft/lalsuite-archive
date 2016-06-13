library("lattice")

zfilter<-function(z, f, ...) {
	x<-filter(Re(z), Re(f), ...)-filter(Im(z), Im(f), ...)
	y<-filter(Re(z), Im(f), ...)+filter(Im(z), Re(f), ...)
	return(x+1i*y)
	}

# Bessel playground

N<-2e6

df<-0.125
df<-0.06

idx<- (1:N)-1

x<-rnorm(N)+1i*rnorm(N)+exp(20*pi*1i*idx/N)

coeffs<-c(-5.32583-0.772804i, -1.31682-0.102469i, -0.583185-0.0306731i)

a<- coeffs[1]*exp(2*pi*1i*idx/N)+coeffs[2]*exp(2*2*pi*1i*idx/N)+coeffs[3]*exp(3*2*pi*1i*idx/N)
mult<-exp(1i*df*(a+Conj(a)))

x2<-x/mult

z<-fft(x)
z2<-fft(x2)


make_bessel_filter3<-function(df, coeffs) {
	if(df<0) {
		df<- -df
		coeffs<- -coeffs
		}

	cnorm<-abs(coeffs)
	cunit<-coeffs/cnorm

	# idx2<- -3:3
	# zcoeffs<-besselJ(cnorm[1]*0.25, idx2)*( (1i*cunit[1])^idx2)
	# zcoeffs[4]<-zcoeffs[4]+besselJ(cnorm[2]*0.25, 0)+besselJ(cnorm[3]*0.25, 0)
	# zcoeffs[2]<-zcoeffs[2]+besselJ(cnorm[2]*0.25, -1)* -1i/cunit[2]
	# zcoeffs[6]<-zcoeffs[6]+besselJ(cnorm[2]*0.25, 1)*1i*cunit[2]
	# zcoeffs[1]<-zcoeffs[1]+besselJ(cnorm[3]*0.25, -1)* -1i/cunit[3]
	# zcoeffs[7]<-zcoeffs[7]+besselJ(cnorm[3]*0.25, 1)*1i*cunit[3]

	idx2<- -3:3
	idx3<- -1:1
	zcoeffs1<-besselJ(cnorm[1]*df*2, idx2)*( (1i*cunit[1])^idx2)
	zcoeffs2<-besselJ(cnorm[2]*df*2, idx3)*( (1i*cunit[2])^idx3)
	zcoeffs3<-besselJ(cnorm[3]*df*2, idx3)*( (1i*cunit[3])^idx3)

	zcoeffs<-zcoeffs1*zcoeffs2[2]
	zcoeffs[1:5]<-zcoeffs[1:5]+zcoeffs1[3:7]*zcoeffs2[1]
	zcoeffs[3:7]<-zcoeffs[3:7]+zcoeffs1[1:5]*zcoeffs2[3]

	zcoeffs4<-zcoeffs

	zcoeffs<-zcoeffs4*zcoeffs3[2]
	zcoeffs[1:4]<-zcoeffs[1:4]+zcoeffs4[4:7]*zcoeffs3[1]
	zcoeffs[4:7]<-zcoeffs[4:7]+zcoeffs4[1:4]*zcoeffs3[3]
	return(zcoeffs)
	}

fold_filter<-function(zcoeffs1, zcoeffs2, k) {
	offset<-floor(length(zcoeffs2)/2)

	z<-zcoeffs1*zcoeffs2[offset+1]

	for(i in 1:offset) {
		m<-k*i
		if(m>=length(z))next
		idx<-1:(length(z)-m)
		z[idx]<-z[idx]+zcoeffs1[idx+m]*zcoeffs2[offset+1-i]
		z[idx+m]<-z[idx+m]+zcoeffs1[idx]*zcoeffs2[offset+1+i]
		}

	return(z)
	}

make_bessel_filter4<-function(df, coeffs, N=4) {
	if(df<0) {
		df<- -df
		coeffs<- -coeffs
		}

	cnorm<-abs(coeffs)
	cunit<-coeffs/cnorm

	print(1i*cunit)

	# idx2<- -3:3
	# zcoeffs<-besselJ(cnorm[1]*0.25, idx2)*( (1i*cunit[1])^idx2)
	# zcoeffs[4]<-zcoeffs[4]+besselJ(cnorm[2]*0.25, 0)+besselJ(cnorm[3]*0.25, 0)
	# zcoeffs[2]<-zcoeffs[2]+besselJ(cnorm[2]*0.25, -1)* -1i/cunit[2]
	# zcoeffs[6]<-zcoeffs[6]+besselJ(cnorm[2]*0.25, 1)*1i*cunit[2]
	# zcoeffs[1]<-zcoeffs[1]+besselJ(cnorm[3]*0.25, -1)* -1i/cunit[3]
	# zcoeffs[7]<-zcoeffs[7]+besselJ(cnorm[3]*0.25, 1)*1i*cunit[3]

	idx2<- -N:N
	idx3<- -N:N
	zcoeffs1<-besselJ(cnorm[1]*df*2, idx2)*( (1i*cunit[1])^idx2)
	zcoeffs2<-besselJ(cnorm[2]*df*2, idx3)*( (1i*cunit[2])^idx3)
	zcoeffs3<-besselJ(cnorm[3]*df*2, idx3)*( (1i*cunit[3])^idx3)

	print(zcoeffs1)
	print(abs(zcoeffs2))
	print(abs(zcoeffs3))
	zcoeffs<-fold_filter(zcoeffs1, zcoeffs2, 2)
	print(zcoeffs)
	zcoeffs<-fold_filter(zcoeffs, zcoeffs3, 3)
	print(zcoeffs)
	return(zcoeffs)
	}

zcoeffs<-make_bessel_filter4(df, coeffs)

z3<-0
for(i in 1:7) {
	z3<-z3+z2[i:(i+200)]*zcoeffs[8-i]
	}

cat(paste("filter_coeffs_gold[", 0:6, "].re=", Re(zcoeffs[1:7]), ";", sep="", collapse="\n"), "\n")
cat(paste("filter_coeffs_gold[", 0:6, "].im=", Im(zcoeffs[1:7]), ";", sep="", collapse="\n"), "\n")
#plot(abs(zfilter(z2, zcoeffs))[1:40])
