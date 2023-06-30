c-- MEM_EST --------------------------------------------------------------------
c
c  Subroutine to make MEM estimate of the directional spectrum
c  This version reads in normalized four. coeff., so no division by a0
c
c  Input:  Energy and 4 normalized directional fourier coeff.
c          from a buoy or slope array.
c
c  Output: Directional spectrum in 5 degree directional bins (in
c           (ds) the "compass heading" coordinate frame).
c
c            72 directional bands. 1=5 degress....72=360 degrees
c            where the direction is the compass heading from which
c            the wave energy is arriving. e.g. 270=from west, 180=
c            from south etc.
c
c  The ds estimate is set to 0 for any band in which direction is not
c  resolved (dir = -1) or where the coeffs could not be normalized.
c-------------------------------------------------------------------------------
      subroutine mem_est(sp_data, mem_out)

        real chk, d(360)
        type(sp_data_block) sp_data
        type(mem_data) mem_out

c-- Initialize number of dir, freq bands

        mem_out%freq_bands = sp_data%bands
        mem_out%dir_bands = MEM_dir_bands

c-- Loop thru freq bands
c-- moments are first four fourier coeff. of the directional
c-- distribution, MEM uses normalized values

        do i= 1, sp_data%bands

          mem_out%freq(i) = sp_data%freq(i)
          mem_out%band_width(i) = sp_data%band_width(i)

c-- get MEM estimate (d) with 1 deg resolution, d is returned in
c-- compass coordinates

          if (sp_data%dir(i) .ne. -1) call mem(sp_data%a1(i), sp_data%a2(i),
     *      sp_data%b1(i), sp_data%b2(i), d, chk)


c-- merge into 5 deg directional bins, multiply by 0.2 to get
c-- units of m^2/Hz-deg

          do j = 5, 355, 5
            if (sp_data%dir(i) .ne. -1) then
              mem_out%ds(j/5,i) =
     *          .2*sp_data%ener_dens(i)*(d(j-2)+d(j-1)+d(j)+d(j+1)+d(j+2))
            else
              mem_out%ds(j/5,i) = 0
            end if

          end do

          if (sp_data%dir(i) .ne. -1) then
            mem_out%ds(72,i) =
     *        0.2*sp_data%ener_dens(i)*(d(358)+d(359)+d(360)+d(1)+d(2))
          else
            mem_out%ds(72,i) = 0
          end if

        end do

        do k = 5, 360, 5
          mem_out%dir(k/5) = k
        end do

        return
      end subroutine


c-- MEM ------------------------------------------------------------------------
c  Maximum entropy method (MEM) for the estimation of a directional
c  distribution from pitch-and-roll data. By default, covers the
c  full 360 degrees with 1.0-degree resolution. (Call MEMD directly for
c  settings other than the default.)
c
c  (Lygre and Krogstad, JPO v16 1986: NOTE - there is a typo in the
c  paper...BOTH exponetials in eq. 13 should have negative signs.
c  This is fixed in Krogstad, 1991- THH has the book)
c
c  variables:
c    a1,b1,b1,b2 = fourier coef. of dir. dist.
c      (normalized ala Long [1980]) - TRIG COORDINATES
c
c    s = MEM estimate of directional spectrum, 1 deg.
c      resolution - COMPASS COORDINATES
c
c    chk = check factor: MEM likes to make narrow spikes for
c      directional spectra if it can (it fits the a1,a2,
c      b1 abd b2's exactly).  If it does this, then the
c      .1 deg stepsize for making the 1 deg resolution
c      estimate is too coarse.  The check factor should
c      be close to 1.
c-------------------------------------------------------------------------------
        subroutine mem(a1,a2,b1,b2,s,chk)
          real    a1, a2, b1, b2, chk, s(360)
          call memd(a1,a2,b1,b2,1,360,1.0,s,chk)
        end subroutine


c-- MEMD -----------------------------------------------------------------------
c  Does the calcs for MEM above. Call directly if do not need full 360 degrees
c  or 1.0-degree resolution. Only returns values for integer degrees; res
c  must be <= 1.0, and the start direction, begin, must be > 0.
c
c  Modified 6/2011, added normalization to handle narrow peaks in hindcast data
c    When running the full directional span (1 to 360), normalized by tot
c    instead of 360.
c-------------------------------------------------------------------------------
        subroutine memd(a1,a2,b1,b2,begin,ndeg,res,s,chk)
          integer begin, dir, ndeg
          real    a1, a2, b1, b2, chk, s(360), offset, res, rn
          complex c1, c2, p1, p2, e1, e2, x, y

          real,parameter:: dr = 0.0174533

          do i = begin, begin+ndeg-1
            if (i .gt. 360) then
              dir = i - 360
            else
              dir = i
            end if
            s(dir)=0.
          end do

c-- switch to Lygre & Krogstad notation

          d1 = a1
          d2 = b1
          d3 = a2
          d4 = b2

          c1=(1.,0)*d1+(0,1.)*d2
          c2=(1.,0)*d3+(0,1.)*d4

          p1=(c1-c2*conjg(c1))/(1.-cabs(c1)**2)
          p2=c2-c1*p1

          x=1.-p1*conjg(c1)-p2*conjg(c2)

c-- sum over 'ndeg' in steps, get distribution with 'res' degree resolution

          tot=0
          offset = 0.5 * (1.0 - res)

          do rn = begin-offset, begin+ndeg-1+offset, res
            a=rn*dr
            e1=(1.,0)*cos(a)-(0,1.)*sin(a)
            e2=(1.,0)*cos(2*a)-(0,1.)*sin(2*a)
            y=cabs((1.,0)-p1*e1-p2*e2)**2

c-- put in proper 1 deg directional band

            ndir=NINT(rn)
            if (ndir .gt. 360) ndir=ndir-360

c-- normalize by 360/(step size) if not running full 360 degrees

            if (ndeg .ne. 360) then
              s(ndir)=s(ndir)+cabs(x/y)/(360./res)
            else
              s(ndir)=s(ndir)+cabs(x/y)
            end if
            tot=tot+cabs(x/y)
          end do

c--  normalize spectrum for full 360 degree run

          if (ndeg .eq. 360) then
            do i = 1, 360
              s(i) = s(i)/tot
            end do
            chk = 1
          else

c-- tot should = 360.  If directional peak is extremely narrow then
c-- 1 deg resolution may be insufficient and tot .ne. 360

            chk=tot/(360./res)
          end if

          return
        end subroutine
