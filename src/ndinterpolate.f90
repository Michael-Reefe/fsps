
! from:
! https://degenerateconic.com/multidimensional-linear-interpolation.html

function ndinterpolate(narg,arg,nent,ent,table)

use sps_vars
implicit none

real(SP) :: ndinterpolate  
integer,intent(in) :: narg  
integer,dimension(narg),intent(in) :: nent  
real(SP),dimension(narg),intent(in) :: arg  
real(SP),dimension(:),intent(in) :: ent  
real(SP),dimension(:),intent(in) :: table

integer,dimension(2**narg) :: index  
real(SP),dimension(2**narg) :: weight  
logical :: mflag,rflag  
real(SP) :: eta,h,x  
integer :: i,ishift,istep,k,knots,lgfile,&  
           lmax,lmin,loca,locb,locc,ndim

!some error checks:  
if (size(ent)/=sum(nent)) &  
    error stop 'size of ent is incorrect.'  
if (size(table)/=product(nent)) &  
    error stop 'size of table is incorrect.'  
if (narg<1) &  
    error stop 'invalid value of narg.'

lmax = 0  
istep = 1  
knots = 1  
index(1) = 1  
weight(1) = 1.0_SP

main: do i = 1 , narg  
    x = arg(i)  
    ndim = nent(i)  
    loca = lmax  
    lmin = lmax + 1  
    lmax = lmax + ndim  
    if ( ndim>2 ) then  
        locb = lmax + 1  
        do  
            locc = (loca+locb)/2  
            if ( x<ent(locc) ) then  
                locb = locc  
            else if ( x==ent(locc) ) then  
                ishift = (locc-lmin)*istep  
                do k = 1 , knots  
                    index(k) = index(k) + ishift  
                end do  
                istep = istep*ndim  
                cycle main  
                else  
                loca = locc  
            end if  
            if ( locb-loca<=1 ) exit  
        end do  
        loca = min(max(loca,lmin),lmax-1)  
        ishift = (loca-lmin)*istep  
        eta = (x-ent(loca))/(ent(loca+1)-ent(loca))  
    else  
        if ( ndim==1 ) cycle main  
        h = x - ent(lmin)  
        if ( h==0.0_SP ) then  
            istep = istep*ndim  
            cycle main  
        end if  
        ishift = istep  
        if ( x-ent(lmin+1)==0.0_SP ) then  
            do k = 1 , knots  
                index(k) = index(k) + ishift  
            end do  
            istep = istep*ndim  
            cycle main  
        end if  
        ishift = 0  
        eta = h/(ent(lmin+1)-ent(lmin))  
    end if  
    do k = 1 , knots  
        index(k) = index(k) + ishift  
        index(k+knots) = index(k) + istep  
        weight(k+knots) = weight(k)*eta  
        weight(k) = weight(k) - weight(k+knots)  
    end do  
    knots = 2*knots  
    istep = istep*ndim  
end do main

ndinterpolate = 0.0_SP  
do k = 1 , knots  
    i = index(k)  
    ndinterpolate = ndinterpolate + weight(k)*table(i)  
end do

end function ndinterpolate 
