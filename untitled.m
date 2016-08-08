clear



%% --
% Vector

		comp_i = (0d0,1d0)
		
		zr = k*w0*w0/2d0
		
		vecTT = duration / sqrt(8d0*log(2d0))
 		
		vect = t+comp_i*zr
		vecz = z+comp_i*zr
		
		vecRsq = x*x + y*y + (z + comp_i*zr)*(z + comp_i*zr)
		vecR = sqrt(vecRsq)
				
		if (aimag(vecR) < 0) then
			vecR = -vecR
		end if
				
		vectau = vect - vecR		
		
		vecf = (1d0 + (comp_i*vectau)/(k*vecTT*vecTT))**2d0-(1d0/(k*k*vecRsq))*(1d0-vect*vecR/(vecTT*vecTT)+comp_i*k*vecR)
        vecg = -vecf + (2d0/(k*k*vecRsq))*(1d0 - vectau*vecR/(vecTT*vecTT) + comp_i*vecR)
        vech = vecf + 1d0/(k*k*vecRsq)
        
        vecB = 1d0 -1d0/(k*zr)+ 1d0/(k*vecTT)**2d0
        vecA0 = sqrt(1d0/(vecB+1d0/(k*zr*k*zr)))
		
		vecp0const = zr*vecA0*E0/(k*k)
		
		vecp0 = vecp0const*exp(-vectau**2d0/(2d0*vecTT*vecTT))*exp(comp_i*k*vectau)
		
		vecA = E0*zr*vecA0*vecp0/(vecp0const*vecR)
		
		E1temp = realpart(vecA*(vecf+x*x*vecg/vecRsq))
		E2temp = realpart(vecA*x*y*vecg/vecRsq)
		E3temp = realpart(vecA*x*vecz*vecg/vecRsq)
		
		B1temp = 0d0
		B2temp = realpart(vecA*vecz*vech/vecR)
		B3temp = realpart(-vecA*y*vech/vecR)