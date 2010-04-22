module priors

	use utils1

contains

!=======================================================================
! Prior distribution functions: r is a uniform deviate from the unit
!
!

! Uniform[0:1]  ->  Delta[x1]

      function DeltaFunctionPrior(r,x1,x2)

      implicit none

      real*8 r,x1,x2,DeltaFunctionPrior

      DeltaFunctionPrior=x1

      return
      end function DeltaFunctionPrior

!=======================================================================
! Uniform[0:1]  ->  Uniform[x1:x2]

      function UniformPrior(r,x1,x2)

      implicit none

      real*8 r,x1,x2,UniformPrior

      UniformPrior=x1+r*(x2-x1)

      return
      end function UniformPrior

!=======================================================================
! Uniform[0:1]  ->  LogUniform[x1:x2]

      function LogPrior(r,x1,x2)

      implicit none

      real*8 r,x1,x2,LogPrior
      real*8 lx1,lx2

      if (r.le.0.0d0) then
       LogPrior=-1.0d32
      else
       lx1=dlog10(x1)
       lx2=dlog10(x2)
       LogPrior=10.d0**(lx1+r*(lx2-lx1))
      endif

      return
      end function LogPrior

!=======================================================================
! Uniform[0:1]  ->  Sin[x1:x2]  (angles in degrees):

      function SinPrior(r,x1,x2)

      implicit none

      real*8 r,x1,x2,SinPrior
      real cx1,cx2,deg2rad
      parameter(deg2rad=0.017453292)

      cx1=cos(x1*deg2rad)
      cx2=cos(x2*deg2rad)
      SinPrior=1.d0*acos(cx1+r*(cx2-cx1))

      return
      end function SinPrior

!=======================================================================
! Uniform[0:1]  ->  Cauchy[mean=x0,FWHM=2*gamma]

      function CauchyPrior(r,x0,gamma)

      implicit none

      real*8 r,x0,gamma,CauchyPrior
      real Pi
      parameter(Pi=3.141592654)

      CauchyPrior=x0+gamma*tan(Pi*(r-0.5))

      return
      end function CauchyPrior

!======================================================================= 
end module priors
