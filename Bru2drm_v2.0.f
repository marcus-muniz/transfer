      PROGRAM BRU2DRM_v2
c
c *********************************************************************
c               PROGRAM BRU2DRM
c
c     Last modification: 1.0 - Noll (10-Jun-2006)
c                          - This program is derived from the BRU2DE_V5.9 with the
c                            modifications implemented in BRU2D_V5.11, except for the
c                            mesh alone motion.
c                          - This variation moves the mesh rigidly with the body.
c
c     Last modification: 2.0 - Noll (20-Jun-2006)
c                          - This version allows the excitation of the displacement
c                            and velocity components of the problem separetely
c
c               PROGRAM BRU2D (Edge Based Data)
c
c     Last modification: 5.9 - Noll (21-Mar-2006)
c                          - The start subroutine has been slightly altered in order
c                            to implement more accurate restart cases.
c                          - The aerodynamic coefficients history are written
c                            in a incremental way.
c
c     Last modification: 5.8 - Noll (16-Mar-2006)
c                          - Now I am convinced that Fort.10 was not so irrelevant and
c                            I am echoing the input in it again. So the new function
c                            assumed by Fort.10 was transfered to Fort.23.
c                          - I also gave the nsave variable part of its old function
c                            back. So, it determines the saving frequency of the intermediary
c                            solutions for restart cases. I also created two new
c                            variables: nintsave, that holds for the saving frequency of
c                            intermediary solutions for plotting; and nhist, which
c                            determines the saving frequency for the aerodynamic 
c                            coefficients history.
c
c     Last modification: 5.7 - Noll (06-Mar-2006)
c                          - The program writes in Fort.22 the CPU time elapsed during
c                            computation. It dows it twice: one after the preprocessing
c                            phase and again when the code is stoped.
c
c     Last modification: 5.6 - Noll (02-Feb-2006)
c                          - I modified the input data file Fort.1 layout in order to
c                            make it easier to modify. I had to change the indat
c                            subroutine also.
c                          - This modification makes Fort.10 irrelevant and I excluded
c                            it from the code.
c                          - In this modified Fort.1 file I included the option of
c                            saving the entire solution for different time steps in
c                            a new Fort.10 file. This is done in stopit subroutine.
c                          - I also decided to follow the piece of advice L.C. Oliveira
c                            left almost 13 years ago. As I believe an implicit solver
c                            is not going to be implemented very soon, I took the call
c                            on boundary subroutine out of the update subroutine. But
c                            I kept the update subroutine once it already controls all
c                            the mesh motion parameters. With this modification, I can
c                            move the call on update after the convergence checks.
c
c     Last modification: 5.5 - Noll (02-Feb-2006)saving frequency for
c                          - I modified the update subroutine so that it does not
c                            perform mesh motion calculation when the motion stops
c                            for the exponentially-shaped pulse.
c                          - I also included the Fort.21 file with mesh deformation
c                            information.
c
c     Last modification: 5.4 - Noll (01-Sep-2005)
c                          - creation of imove = 2 and 3 options corresponding to the step
c                          and unit sample excitations.
c     
c     Last modification: 5.3 - Noll (15-Jun-2004)
c                          - I removed the nsave2 definition to use nsave
c                            directly
c
c     Last modification: 5.2 - Noll (07-Jun-2004)
c                          - I modified the preproc subroutine so it
c                            considerates the flat plate case
c
c     Last modification: 5.1 - Noll (25-May-2004)
c                         - output the steady Cp distribution on the wall
c                           in Fort.20 and the unsteady to Fort.19 instead
c                           of doing it in Fort.10.                        
c                         - output these results in an ordered fashion to
c                           ease the plotting using Tecplot.
c                         - creation of subroutine pressure
c
c     Last modification: 5.0 - Fred (06-Oct-2003)
c                        - Implementacao da geometric conservation law e de
c                          calculo de velocidade dos nos por segunda ordem.
c                          Coloquei novas variaveis no common e mudei as rotinas
c                          Dynmesh, Meshio, Indat e Update. Criei nova subrotina
c                          GCL.
c                        - A dissipacao do Mavriplis foi mudada para que o modulo
c                          fique dentro do somatorio.
c                        - WARNING: restart nao-estacionario nao foi refeito para
c                          adequar com novas mudancas.
c
c     Last modification: 4.0 - Fred (17/Sep/2003)
c                        - Alguns commons internos que estavam comentados
c                          foram apagados
c                        - Mudanca na rotina de condicao de contorno na parede.
c                          A rotina faz essencialmente a mesma coisa, porem
c                          com menor numero de operacoes
c                        - Dissipacao artificial do Mavriplis refeita do zero,
c                          variavel Abig retirada do common.
c                        - Output de Alpha instantaneo para tela
c                        - Rotina Stopit agora tambem escreve entropia em
c                          torno do perfil, para caso estacionario.
c
c     Last modification: 3.0 - Fred (02/Sep/2003)
c                        - Combinacao das versoes 1.3, 1.4 e 2.
c
c     Last modification: 1.4 - Fred (20/Apr/1999)
c                        - change in the artificial dissipation to
c                        Mavriplis ideas
c
c     Last modification: 1.3 - Fred (29/Mar/1999)
c                        - another new aerodynamic coefficients
c                        calculation routine
c                        - without alpha output to screen
c     
c     Last modification: 1.2 - Fred (27/Mar/1999)
c                        - new aerodynamic coefficients calculation
c                        routine
c                        - instantaneous CP output to compare against
c                        experimental data
c                        - alpha output to screen
c
c     Last modification: 1.1 - Fred (19/Oct/1998)
c                        - files Fort.n.dat moved to just Fort.n
c
c     Last modification: 2.0 - Fred (16/Mar/1998)
c                        - RK 4-stages
c
c     Last modification: 1.0 - Fred (06/Aug/1997)
c                        - unsteady version
c
c     Last modification: Joao Azevedo     Feb./02/96
c
c     This is version 4.1 of the code.
c
c *********************************************************************
c
c     File 'Common_rm_v2.0.f' contains PARAMETER, IMPLICIT and COMMON definitions
c
      include 'Common_rm_v2.0.f'
c
c     Open all the necessary input/output files (optional files are
c     opened at the appropriate subroutines).
c
c     1  -> Input data
c     2  -> Grid file
c     3  -> Connectivity table
c     4  -> Neighbors information
c     5  -> Ghost volumes information
c     7  -> Results of a previous steady-state run
c     8  -> Steady-state conserved variables (Q's)
c     9  -> Maximum residues file
c     10 -> Echo of input data
c     11 -> L2 residues
c     14 -> Aerodynamic Coefficients
c     15 -> Results of a previous time-accurate run
c     16 -> Time-accurate conserved variables
c     19 -> Unsteady Cp's
c     20 -> Steady Cp's
c     21 -> Mesh deformation
c     22 -> CPU time
c     23 -> Solutions saved at a desired frequency
c
      open(1,file='Fort.1')
      open(2,file='Fort.2')
      open(3,file='Fort.3')
      open(4,file='Fort.4')
      open(5,file='Fort.5')
      open(9,file='Fort.9')
      open(10,file='Fort.10')
      open(11,file='Fort.11')
      open(14,file='Fort.14')
c
c     Read input data and set up initial conditions.
c
      call indat
c
c     Set up iteration counter to start main loop. For the unsteady
c     flow case, we perform another 'steady' iteration before going
c     to move the mesh. This is enforced into the convenient sub-
c     routines.
c
      niter = nstart
      nmiter = nstart + nmiter
c
c     Start Main Loop.
c
   20 continue
      niter = niter + 1
c
c     Compute current residue in the whole field. This subroutine
c     also returns the information on the maximum residue in the
c     field which is important for convergence test.
c
        print*, "residue"
      call residue
      if (iflow2.eq.0) then
        if(rhsmax.lt.epscon) call stopit(1)
        if(rhsmax.gt.epsblw) call stopit(2)
      endif
c
c     Perform one time step using the 4 or 5-stage Runge-Kutta explicit
c     time-stepping scheme.
c
      if (itime.eq.0) then
         call rk5s
      elseif (itime.eq.1) then
         call rk4s
      endif
c
c     The "ghost" volumes are taken care of by subroutine Boundary.
c
      call boundary
c
c     Calculate the aerodynamic coefficients.
c
      if(mod(niter,nhist).eq.0) then
         call history
      endif
c
c     Save the wall pressure distribution at experimental alphas.
c     This has to be changed so that the user can determine at which
c     alphas he wants it to be saved. But that will be done sometime
c     in the future.
c
c      alphai = alphat*180.d0/3.1415926d0 + alpha
c      if((alphai.gt.1.089  .and. alphai.lt.1.091) .or.
c     &   (alphai.gt.2.339  .and. alphai.lt.2.341) .or.
c    &   (alphai.gt.2.009  .and. alphai.lt.2.011) .or.
c     &   (alphai.gt.0.519  .and. alphai.lt.0.521) .or.
c     &   (alphai.gt.-1.251 .and. alphai.lt.-1.249) .or.
c     &   (alphai.gt.-2.411 .and. alphai.lt.-2.409) .or.
c     &   (alphai.gt.-2.001 .and. alphai.lt.-1.999) .or.
c     &   (alphai.gt.-0.541 .and. alphai.lt.-0.539)) then
c         if(niter.ne.nmiter) call stopit(4)
c      end if
c
c     Check if intermediary solutions should be saved for plotting.
c
      if(iresp .eq. 1) then
         if(mod(niter,nintsave).eq.0) call intermed
      endif
c
c     Check if the pressure distribution is also to be saved
c
      if(ipress .eq. 1) then
         if(mod(niter,npressave).eq.0) call pressure
      end if
c
c     Update interior volume properties, apply boundary conditions
c     and, for the unsteady flow case, move the mesh.
c
      call update
c
c     Save the current solution in the disk every "nsave" iteration.
c
      if(mod(niter,nsave).eq.0) then
         if(niter.ne.nmiter) call stopit(3)
      endif
c
c     If maximum number of iterations has been exceeded, stop the
c     code. Otherwise, continue the iteration.
c
      if(niter.lt.nmiter) go to 20
      call stopit(0)
      stop
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine BOUNDARY
c
c     This subroutine computes the properties of the "ghost", or
c     "slave" control volumes. This is the actual implementation of
c     boundary conditions in the present finite volume context (i.e.,
c     cell-centered).
c     Here, I have specified a few types of boundary conditions. The
c     type of boundary that a given ghost volume belongs to is
c     specified in its share of information in the ghost volume
c     connectivity table. If the information there indicates a type
c     of b.c. not specified, the program will bomb out (with the
c     appropriate message, of course).
c     UnFortunately, I don't think this subroutine will be much
c     vectorizable. In order to make things general, I will pretty
c     much have to do it volume by volume.
c
      subroutine boundary
c
      include 'Common_rm_v2.0.f'
c
      real*8 dltx, dlty, dltxy, dltx2, dlty2, dlts2
     &       u, v, ugh, vgh
c
c     Loop over all the "ghost" volumes.
c
      do 500 ig=1,nvolg
c
c        Get information from connectivity table and check type of
c        boundary.
c
         ibtype = ighost(ig,1)
         intngb = ighost(ig,2)
         in1 = ighost(ig,3)
         in2 = ighost(ig,4)
c
c        Boundary type # 1 is a solid wall boundary.
c
         if(ibtype.eq.1) then
c
c-----------------------------------------------------------------------
c        Este eh o jeito que o Luis Claudio implementou
c        (12 somas, 16 multiplicacoes, 1 raiz, 20 atribuicoes)
c
c            aux1 = xynode(in2,1) - xynode(in1,1)
c            aux2 = xynode(in2,2) - xynode(in1,2)
c            d = dsqrt(aux1*aux1 + aux2*aux2)
c            sx = aux2 / d
c            sy = aux1 / d
c            if(nvisc.eq.0) then
c               u = q(intngb,2) / q(intngb,1)
c               v = q(intngb,3) / q(intngb,1)
c               uu = u-uvol(intngb)
c               vv = v-vvol(intngb)
c               uun = sx*uu - sy*vv
c               uut = sy*uu + sx*vv
c               uun = - uun
c               uu = sx*uun + sy*uut
c               vv = - sy*uun + sx*uut
c               u = uu + uvol(intngb)
c               v = vv + vvol(intngb)
c               rghost(ig,1) = q(intngb,1)
c               rghost(ig,2) = rghost(ig,1) * u
c               rghost(ig,3) = rghost(ig,1) * v
c               rghost(ig,4) = q(intngb,4)
c            else
c               rghost(ig,1) = q(intngb,1)
c               rghost(ig,2) = - q(intngb,2)
c               rghost(ig,3) = - q(intngb,3)
c               rghost(ig,4) = q(intngb,4)
c            end if
c            go to 490
c-----------------------------------------------------------------------
c      Jeito do Fred, somente um pouquinho mais otimizado (08-sep-2003)
c      (9 somas, 15 multiplicacoes, nenhuma raiz, 14 atribuicoes)
c
            dltx = xynode(in2,1)-xynode(in1,1)
            dlty = xynode(in2,2)-xynode(in1,2)
            dltxy = dltx*dlty
            dltx2 = dltx*dltx
            dlty2 = dlty*dlty
            dlts2 = dltx2+dlty2
c
            if(nvisc.eq.0) then
c
c           If nvisc = 0, this is an Euler computation. Here we take
c           into account the unsteady flow case by using the components
c           of the cell velocity.
c
               u = q(intngb,2)/q(intngb,1) - uvolg(ig)
               v = q(intngb,3)/q(intngb,1) - vvolg(ig)
c
               ugh = (2.0d0*v*dltxy + u*(dltx2-dlty2)) / dlts2
               vgh = (2.0d0*u*dltxy + v*(dlty2-dltx2)) / dlts2
c
               rghost(ig,1) = q(intngb,1)
               rghost(ig,2) = rghost(ig,1) * (ugh + uvolg(ig))
               rghost(ig,3) = rghost(ig,1) * (vgh + vvolg(ig))
               rghost(ig,4) = q(intngb,4)
c
            else
c
c           If nvisc = 1 or 2, this is a Navier-Stokes computation
c           (I'm assuming here an adiabatic wall case).
c           Here, as it is, I am certainly getting the correct
c           momentum components at the wall, but I'm not sure about
c           the energy... It seems to me that something fancier must
c           be done with the energy to get it right at wall for the
c           N-S case, but I'll not think about it today.
c
               rghost(ig,1) = q(intngb,1)
               rghost(ig,2) = - q(intngb,2)
               rghost(ig,3) = - q(intngb,3)
               rghost(ig,4) = q(intngb,4)
c
            end if
            go to 490
         end if
c
c        Boundary type # 2 is a freestream boundary.
c
c        This is the old way of implementing this boundary condition,
c        that is, imposing freestream at the far field.
c
c         if(ibtype.eq.2) then
c            rghost(ig,1) = 1.0d0
c            rghost(ig,2) = rghost(ig,1) * fsmach * dcos(alpha)
c            rghost(ig,3) = rghost(ig,1) * fsmach * dsin(alpha)
c            rghost(ig,4) = (1.0d0 / (gamma*(gamma - 1.0d0)))
c     &                   + 0.5d0*rghost(ig,1)*fsmach*fsmach
c            go to 490
c         end if
c
c
c        Here comes the new way, i.e., using Riemann invariants.
c        It is also important to observe that a test for supersonic
c        normal velocity component at the boundary should be
c        implemented in the future. This is actually very easy to
c        do, but I didn't want to use the time to do it now.
c
         if(ibtype.eq.2) then
c
c           Compute direction cosines.
c
            aux1 = xynode(in2,1) - xynode(in1,1)
            aux2 = xynode(in2,2) - xynode(in1,2)
            d = dsqrt(aux1*aux1 + aux2*aux2)
            sx = aux2 / d
            sy = aux1 / d
c
c           Compute incoming and outgoing Riemann invariants.
c
            uinf = fsmach * dcos(alpha)
            vinf = fsmach * dsin(alpha)
            qninf = sx*(uinf - uvolg(ig)) - sy*(vinf - vvolg(ig))
            rrinf = qninf - (2.0d0 / (gamma - 1.0d0))
c
            ue = q(intngb,2) / q(intngb,1)
            ve = q(intngb,3) / q(intngb,1)
            qne = sx*(ue - uvolg(ig)) - sy*(ve - vvolg(ig))
            pe = (gamma - 1.0d0) * (q(intngb,4) - 0.5d0 *
     &           q(intngb,1) * (ue*ue + ve*ve))
            ae = dsqrt(gamma * pe / q(intngb,1))
            rre = qne + ((2.0d0 / (gamma - 1.0d0)) * ae)
c
c           Compute normal velocity component and speed of sound
c           at the interface.
c
            qn = 0.5d0 * (rre + rrinf)
            a = ((gamma - 1.0d0) / 4.0d0) * (rre - rrinf)
c
c           Compute the other properties at the interface.
c
            if(qn.gt.0) then
c
c              This is the outflow case.
c
               aux1 = 1.0d0 / (gamma - 1.0d0)
               rho = ((q(intngb,1)**gamma * a*a)/(gamma * pe))**aux1
               u = ue + (qn - qne) * sx
               v = ve - (qn - qne) * sy
c
c              This is the inflow case.
c
            else
               aux1 = 2.0d0 / (gamma - 1.0d0)
               rho = a**aux1
               u = uinf + (qn - qninf) * sx
               v = vinf - (qn - qninf) * sy
            end if
c
            p = (rho * a * a) / gamma
            e = (p / (gamma - 1.0d0)) + (0.5d0 * rho * (u*u + v*v))
c
c           Finally, compute the ghost volume information.
c
            rghost(ig,1) = (2.0d0 * rho) - q(intngb,1)
            rghost(ig,2) = (2.0d0 * rho * u) - q(intngb,2)
            rghost(ig,3) = (2.0d0 * rho * v) - q(intngb,3)
            rghost(ig,4) = (2.0d0 * e) - q(intngb,4)
c
            go to 490
         end if
c
c        Boundary type # 3 is an entrance boundary.
c
         if(ibtype.eq.3) then
c
c           Determine whether flow is subsonic or supersonic.
c
            u = q(intngb,2) / q(intngb,1)
            v = q(intngb,3) / q(intngb,1)
            p = (gamma - 1.0d0) * (q(intngb,4) - 0.5d0*q(intngb,1)*
     &                                                 (u*u + v*v))
            a = dsqrt((gamma*p)/q(intngb,1))
c
c           For subsonic flow, the cartesian u component is
c           extrapolated, and the rest must be given.
c
            if((u-uvolg(ig)).lt.a) then
               v = u * dtan(alpha)
               aux1 = 1.0d0 - (((gamma - 1.0d0)/(gamma + 1.0d0))*
     &                         (u*u + v*v))
               eint = ((gamma + 1.0d0)/(2.0d0*gamma*(gamma - 1.0d0)))
     &                * aux1
               p = ((gamma + 1.0d0)/(2.0d0*gamma))*(aux1**(gamma/
     &             (gamma - 1.0d0)))
               rghost(ig,1) = p / ((gamma - 1.0d0)*eint)
               rghost(ig,2) = rghost(ig,1) * u
               rghost(ig,3) = rghost(ig,1) * v
               rghost(ig,4) = rghost(ig,1) * (eint + 0.5d0*(u*u + v*v))
c
c           For supersonic flow, all quantities at the entrance
c           must be given.
c
            else
               rghost(ig,1) = 1.0d0
               rghost(ig,2) = rghost(ig,1) * fsmach * dcos(alpha)
               rghost(ig,3) = rghost(ig,1) * fsmach * dsin(alpha)
               rghost(ig,4) = (1.0d0 / (gamma*(gamma - 1.0d0)))
     &                      + 0.5d0*rghost(ig,1)*fsmach*fsmach
            end if
            go to 490
         end if
c
c        Boundary type # 4 is an exit boundary.
c
         if(ibtype.eq.4) then
c
c           Determine whether flow is subsonic or supersonic.
c
            u = q(intngb,2) / q(intngb,1)
            v = q(intngb,3) / q(intngb,1)
            p = (gamma - 1.0d0) * (q(intngb,4) - 0.5d0*q(intngb,1)*
     &                                                 (u*u + v*v))
            a = dsqrt((gamma*p)/q(intngb,1))
c
c           For subsonic flow, the exit pressure is fixed and all
c           other quantities are obtained from interior information.
c           I'm assuming that the pressure is fixed as a fraction of
c           the appropriate freestream static pressure or "reservoir"
c           stagnation pressure. Hence, this condition will also
c           depend on the "flow type".
c
            if((u-uvolg(ig)).lt.a) then
c
               if(iflow.eq.3) then
                  p = pexit*((gamma + 1.0d0)/(2.0d0*gamma))
               else
                  p = pexit / gamma
               end if
c
               rghost(ig,1) = q(intngb,1)
               rghost(ig,2) = q(intngb,2)
               rghost(ig,3) = q(intngb,3)
               rghost(ig,4) = p/(gamma - 1.0d0) + 0.5d0*(rghost(ig,2)*
     &          rghost(ig,2) + rghost(ig,3)*rghost(ig,3))/rghost(ig,1)
c
c           For supersonic exit, all boundary quantities are obtained
c           by extrapolation of interior information.
c
            else
               rghost(ig,1) = q(intngb,1)
               rghost(ig,2) = q(intngb,2)
               rghost(ig,3) = q(intngb,3)
               rghost(ig,4) = q(intngb,4)
            end if
            go to 490
         end if
c
c        Boundary type # 5 is a symmetry boundary.
c
         if(ibtype.eq.5) then
            aux1 = xynode(in2,1) - xynode(in1,1)
            aux2 = xynode(in2,2) - xynode(in1,2)
            d = dsqrt(aux1*aux1 + aux2*aux2)
            sx = aux2 / d
            sy = aux1 / d
c
            u = q(intngb,2) / q(intngb,1)
            v = q(intngb,3) / q(intngb,1)
            un = sx*u - sy*v
            ut = sy*u + sx*v
            un = - un
            u = sx*un + sy*ut
            v = - sy*un + sx*ut
c
            rghost(ig,1) = q(intngb,1)
            rghost(ig,2) = rghost(ig,1) * u
            rghost(ig,3) = rghost(ig,1) * v
            rghost(ig,4) = q(intngb,4)
            go to 490
         end if
c
c        If the code ever gets to this point, something is wrong.
c        This means that the type of boundary specified for this
c        particular volume is not one of the allowed possibilities.
c        Hence, print some error message and stop the code.
c
c         write(*,901) ig,ibtype
         stop
c
  490    continue
  500 continue
c
      return
c
c     Formats.
c
  901 format(///,1x,'BOUNDARY, ghost volume = ',i6,', boundary type',
     & ' specified is not allowed.',//,9x,'(ibtype = ',i6,')')
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine CONVEC
c
c     This subroutine computes the convective operator, C(Q), for
c     all "real" control volumes.
c
c     It could happen that this subroutine will be rather difficult
c     to vectorize, because I'll be doing quite a bit of testing and
c     searching in the connectivity table. I'll try to do my best,
c     though !! (This is an old comment. I am not sure it is still
c     valid with the new way of doing things, but ...)
c
      subroutine convec
c
      include 'Common_rm_v2.0.f'
c
c     Zero out the convective operator array.
c
      do 30 n=1,4
         do 25 i=1,nvol
            cpen(i,n) = 0.d0
   25    continue
   30 continue
c
c     Since the code can actually run different spatial discretization
c     schemes, first of all we must decide which scheme the user has
c     selected.
c
      if(ischeme.gt.1) go to 300
c
c     If we get to this point, it means that we are running the standard
c     Jameson algorithm (central difference type algorithm). I emphasize
c     that here I am implementing the original idea of a cell centered,
c     central difference-type scheme, that is, I compute interface
c     averages of the conserved variables and then form the fluxes.
c
c     In some sense, all that one has to do is to loop over all edges
c     and accumulate the contributions to the convective operator.
c
      do 250 j=1,nsides
c
c        Get the information from the edge-based database array.
c
         in1 = iedges(j,1)
         in2 = iedges(j,2)
         i = iedges(j,3)
         nb = iedges(j,4)
c
c        Decide whether the neighbor is a "real" volume or a "ghost"
c        volume.
c
         if(nb.lt.0) go to 150
c
c        For "real" volumes, compute the average of the conserved
c        variables.
c
         do 100 n=1,4
            qav(n) = 0.5d0 * (q(i,n) + q(nb,n))
  100    continue
c
c        Compute the edge velocity by means of an average of nodes
c        velocities.
c
         uedge = vside(j,1)
         vedge = vside(j,2)
c
c        Compute the E and F fluxes at the interface.
c
         call flux
c
c        Now, add the contribution of this edge to the flux balance
c        of the i-th control volume and subtract it from the flux
c        balance of the nb-th control volume (remember, this is the
c        case in which nb is necessarily an interior volume).
c
         do 120 n=1,4
            aux = evec(n) * (xynode(in2,2) - xynode(in1,2))
     &          - fvec(n) * (xynode(in2,1) - xynode(in1,1))
            cpen(i,n) = cpen(i,n) + aux
            cpen(nb,n) = cpen(nb,n) - aux
  120    continue
         go to 240
c
c        For "ghost" volumes, you just have to do the same only
c        remembering that the neighbor now is in the ghost volume
c        list. Its location on the ghost volume table is (-nb).
c        The exceptions to this rule occur when the boundary is a
c        solid wall boundary or a symmetry boundary. In these cases
c        we have to enforce the zero convective flux condition. I must
c        add that this may be a bit of over-carefulness, because the
c        b.c.'s should be already taking care of it. But, what the
c        heck ??! Better safe than sorry !!
c
  150    continue
         nb = - nb
         ibtype = ighost(nb,1)
         if(ibtype.eq.1.or.ibtype.eq.5) go to 220
c
c        For all other boundaries, just do the same as with the
c        interior volumes.
c
c        Compute the average of the conserved variables.
c
         do 170 n=1,4
            qav(n) = 0.5d0 * (q(i,n) + rghost(nb,n))
  170    continue
c
c        Compute the edge velocity by means of an average of nodes
c        velocities.
c
         uedge = vside(j,1)
         vedge = vside(j,2)
c
c        Compute the E and F fluxes at the interface.
c
         call flux
c
c        Now, add the contribution of this edge to the flux balance
c        of the i-th control volume. There is no flux balance equation
c        for the ghost volume.
c
         do 200 n=1,4
            aux = evec(n) * (xynode(in2,2) - xynode(in1,2))
     &          - fvec(n) * (xynode(in2,1) - xynode(in1,1))
            cpen(i,n) = cpen(i,n) + aux
  200    continue
         go to 240
c
c        For wall and symmetry boundaries, zero the convective flux
c        contribution.
c
  220    continue
c
c        Compute the average of the conserved variables.
c
         do 230 n=1,4
            qav(n) = 0.5d0 * (q(i,n) + rghost(nb,n))
  230    continue
c
c        Compute the edge velocity by means of an average of nodes
c        velocities.
c
         uedge = vside(j,1)
         vedge = vside(j,2)
c
c        Then, compute the interface pressure and the flux.
c
         const = gamma - 1.0d0
         pint = const * (qav(4) - 0.5d0 * (qav(2)*qav(2) +
     &                                     qav(3)*qav(3)) / qav(1))
c
         cpen(i,2) = cpen(i,2)+pint*(xynode(in2,2) - xynode(in1,2))
         cpen(i,3) = cpen(i,3)-pint*(xynode(in2,1) - xynode(in1,1))
         cpen(i,4) = cpen(i,4) +
     &               pint*uedge*(xynode(in2,2) - xynode(in1,2)) -
     &               pint*vedge*(xynode(in2,1) - xynode(in1,1))
c
  240    continue
  250 continue
c
      return
c
c     WARNING: Only Jameson centered scheme will work for
c               unsteady cases ...
c
c     Here we are doing (or sort of) the other schemes. The possible
c     options are: (2) MacCormack's scheme, (3) 1st-order Van Leer
c     flux-vector splitting (FVS) scheme, (4) 1st-order accurate
c     Liou-Steffen FVS algorithm, (5) 2nd-order Van Leer FVS scheme,
c     and (6) 2nd-order Liou-Steffen FVS scheme.
c
c     Options (2), (4) and (6) have not been implemented yet.
c
  300 continue
c
c     MacCormack's scheme is actually done somewhere else. Therefore,
c     if we get to this point with ischeme=2, there is a problem.
c     Hence, print error message and stop.
c
      if(ischeme.eq.2) then
         write(*,901)
         write(*,951)
         stop
      end if
c
c     Implementation of the 1st-order Van Leer FVS scheme.
c
      if(ischeme.eq.3) then
c
c        First, one has to generate the "plus" and "minus" fluxes for
c        all control volumes, including the ghost ones.
c        This will be done in a separated routine just to make things
c        neater. I emphasize that this is a calculation over control
c        volumes, not interfaces (or edges).
c
c
c        Probably, the Liou-Steffen case can be implemented together
c        with the Van Leer scheme. Simply check which is the scheme
c        selected here, and call the appropriate flux calculation
c        routine, i.e., VLflux or LSflux.
c
c
         call vlflux
c
c        Now, just loop over the edges, compute the interface flux
c        contributions and accumulate them on the appropriate
c        convective operators.
c
         do 550 j=1,nsides
c
c           Get the information from the edge-based database array.
c
            in1 = iedges(j,1)
            in2 = iedges(j,2)
            i = iedges(j,3)
            nb = iedges(j,4)
c
c           Decide whether the neighbor is a "real" volume or a
c           "ghost" volume.
c
            if(nb.lt.0) go to 450
c
c           For "real" volumes, simply compute the interface flux in
c           a straightforward manner. Just remember that:
c           F_interface = F_plus(left) + F_minus(right).
c           However, there is the problem of deciding who is "left"
c           and who is "right".
c           The interface flux contribution is added to the balance
c           of the i-th control volume and subtracted from the flux
c           balance of the nb-th control volume (remember, this is the
c           case in which nb is necessarily an interior volume).
c
            dxik = xynode(in2,1) - xynode(in1,1)
            dyik = xynode(in2,2) - xynode(in1,2)
            do 420 n=1,4
c
c              Compute the interface E flux.
c
               if(dyik.ge.0.0d0) then
                  aux1 = eplus(i,n) + eminus(nb,n)
               else
                  aux1 = eplus(nb,n) + eminus(i,n)
               end if
c
c              Compute the interface F flux.
c
               if(dxik.le.0.0d0) then
                  aux2 = fplus(i,n) + fminus(nb,n)
               else
                  aux2 = fplus(nb,n) + fminus(i,n)
               end if
c
c              Sum them up.
c
               aux = (aux1 * dyik) - (aux2 * dxik)
               cpen(i,n) = cpen(i,n) + aux
               cpen(nb,n) = cpen(nb,n) - aux
  420       continue
            go to 540
c
c           For "ghost" volumes, you just have to do the same only
c           remembering that the neighbor now is in the ghost volume
c           list. Its location on the ghost volume table is (-nb).
c
c           I tried in the past to strictly enforce a no convective
c           flux condition at wall and symmetry boundaries.
c           The implementation used for those boundaries was exactly
c           similar to what I do in the centered scheme case.
c           But, this did not work. So, I am treating all boundaries
c           in the same fashion, that is, exactly as I treat interior
c           interfaces.
c
  450       continue
            nb = - nb
c
c           For all boundaries, just do the same as with the interior
c           volumes.
c           Add the contribution of this edge to the flux balance of
c           the i-th control volume. There is no flux balance equation
c           for the ghost volume.
c
            dxik = xynode(in2,1) - xynode(in1,1)
            dyik = xynode(in2,2) - xynode(in1,2)
            do 500 n=1,4
c
c              Compute the interface E flux.
c
               if(dyik.ge.0.0d0) then
                  aux1 = eplus(i,n) + egminus(nb,n)
               else
                  aux1 = egplus(nb,n) + eminus(i,n)
               end if
c
c              Compute the interface F flux.
c
               if(dxik.le.0.0d0) then
                  aux2 = fplus(i,n) + fgminus(nb,n)
               else
                  aux2 = fgplus(nb,n) + fminus(i,n)
               end if
c
c              Sum them up.
c
               aux = (aux1 * dyik) - (aux2 * dxik)
               cpen(i,n) = cpen(i,n) + aux
  500       continue
c
c
  540       continue
  550    continue
c
         return
c
      end if
c
c
c     Implementation of the 1st-order Liou-Steffen FVS scheme.
c
      if(ischeme.eq.4) then
c
c        The Liou-Steffen FVS scheme has not been implemented yet.
c        Hence, print error message and stop.
c
         write(*,951) ischeme
         stop
c
      end if
c
c     Implementation of the 2nd-order Van Leer FVS scheme.
c
      if(ischeme.eq.5) then
c
c        The subroutine VLFLUX2 is the one which actually implements
c        the 2nd order Van Leer FVS scheme. It actually does all the
c        necessary calculations and, therefore, the code here is very
c        simple. Only remember that, when running this option, you
c        must make sure that subroutine INDAT calls subr. PREPROC2,
c        because this is the one which does all the geometrical
c        calculations that VLFLUX2 is going to need. Moreover, VLFLUX2
c        already returns the interface fluxes E_ik and F_ik.
c
c        This is going to be a simple loop over edges to accumulate the
c        contributions for the convective operator.
c
         do 750 j=1,nsides
c
c           Get the information from the edge-based database array.
c
            in1 = iedges(j,1)
            in2 = iedges(j,2)
            i = iedges(j,3)
            nb = iedges(j,4)
c
c           Perform some preliminary calculations.
c
            dxik = xynode(in2,1) - xynode(in1,1)
            dyik = xynode(in2,2) - xynode(in1,2)
c
c           VLflux2 does all the calculation for this edge.
c
            call vlflux2(j)
c
c           Decide whether the neighbor is a "real" volume or a "ghost"
c           volume.
c
            if(nb.ge.0) then
c
c              This is the case in which neighbor is a real volume.
c              Therefore, add the contribution of this edge to the flux
c              balance of the i-th control volume and subtract it from
c              the flux balance of the nb-th control volume (remember,
c              this is the case in which nb is necessarily an interior
c              volume).
c
               do 720 n=1,4
                  aux = evec(n) * dyik - fvec(n) * dxik
                  cpen(i,n) = cpen(i,n) + aux
                  cpen(nb,n) = cpen(nb,n) - aux
  720          continue
c
            else
c
c              This is the case in which neighbor is a ghost volume.
c              Therefore, simply add the contribution of this edge to
c              the flux balance of the i-th control volume.
c
               do 740 n=1,4
                  aux = evec(n) * dyik - fvec(n) * dxik
                  cpen(i,n) = cpen(i,n) + aux
  740          continue
c
            end if
c
c           It is that simple !!!  We are done with this edge.
c
  750    continue
c
         return
c
      end if
c
c     The 2nd-order Liou-Steffen FVS scheme has not been implemented
c     yet.
c     Hence, print error message and stop.
c
      write(*,951) ischeme
      stop
c
c     Formats.
c
  901 format(///,5x,'You selected MacCormack scheme and, somehow, ',/,
     & 5x,'the code ended up in the wrong routine (CONVEC).')
  951 format(///,5x,'You have selected an option which has not been ',
     & 'implemented yet.',/,5x,'Valid options are: (1) Jameson, ',
     & '(3) Van Leer FVS, ',/,5x,'and (5) Van Leer FVS - 2nd order',/,
     & 5x,'You have selected option no. ',i6)
c
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine DISSIP
c
c     This subroutine computes the dissipation operator for all the
c     "real" volumes. This subroutine might also have a hard time to
c     get substantial vectorization.
c
c
c     The fact that we now have an edge-based data structure has not
c     been incorporated into this routine yet. This was not done
c     because I did not want to invest the time to do it now. However,
c     parts of this routine (namely, the loop labeled 300) could have
c     their computational time cut in half if the edge-based structure
c     is used.
c
c
      subroutine dissip
c
      include 'Common_rm_v2.0.f'
c
c
c     Define some temporary storage.
c
      dimension pp(nmaxvl),ppg(nmaxgh),anu(nmaxvl),anabla2(nmaxvl,4)
      dimension dik2(4),dik4(4)
c
c     Zero out the artificial dissipation operator array.
c
      do 40 n=1,4
         do 30 i=1,nvol
            dpen(i,n) = 0.d0
   30    continue
   40 continue
c
c     Probably, this would be the place to introduce a test to check
c     which scheme is being used in the current computation.
c     The idea is that I would keep the RK5S routine as it is, because
c     that routine is already fairly complicated with several possible
c     options. Hence, I would simply let the scheme zero out the
c     artificial dissipation operator and add zeros to the flux balance
c     in the upwind cases.
c     I would be paying an extra cost here, but probably it is worth in
c     order to keep the code cleaner for future modifications.
c
      if(ischeme.gt.2) return
c
c     Loop over all real volumes and compute the pressure.
c
      const = gamma - 1.0d0
      do 70 i=1,nvol
         pp(i) = const * (q(i,4) - 0.5d0 * (q(i,2)*q(i,2)
     &                                    + q(i,3)*q(i,3)) / q(i,1))
   70 continue
c
c     Loop over all ghost volumes and compute the pressure.
c
      const = gamma - 1.0d0
      do 80 ig=1,nvolg
         ppg(ig) = const * (rghost(ig,4) - 0.5d0 * (rghost(ig,2)*
     &     rghost(ig,2) + rghost(ig,3)*rghost(ig,3)) / rghost(ig,1))
   80 continue
c
c     Zero out the auxiliary array which will accumulate the undivided
c     Laplacian of the conserved quantities.
c
      do 100 n=1,4
         do 90 i=1,nvol
            anabla2(i,n) = 0.d0
   90    continue
  100 continue
c
c     Loop over all real volumes and compute the pressure switch term
c     and the undivided Laplacian of the conserved variables.
c
      do 200 i=1,nvol
c
c        Get neighborhood information.
c
         nb1 = neighbor(i,1)
         nb2 = neighbor(i,2)
         nb3 = neighbor(i,3)
c
c        Get information from 1st neighbor.
c
         if(nb1.gt.0) then
            pp1 = pp(nb1)
            do 125 n=1,4
               anabla2(i,n) = anabla2(i,n) + q(nb1,n)
  125       continue
         else
            nb1 = - nb1
            pp1 = ppg(nb1)
            do 130 n=1,4
               anabla2(i,n) = anabla2(i,n) + rghost(nb1,n)
  130       continue
         end if
c
c        Get information from 2nd neighbor.
c
         if(nb2.gt.0) then
            pp2 = pp(nb2)
            do 135 n=1,4
               anabla2(i,n) = anabla2(i,n) + q(nb2,n)
  135       continue
         else
            nb2 = - nb2
            pp2 = ppg(nb2)
            do 140 n=1,4
               anabla2(i,n) = anabla2(i,n) + rghost(nb2,n)
  140       continue
         end if
c
c        Get information from 3rd neighbor.
c
         if(nb3.gt.0) then
            pp3 = pp(nb3)
            do 145 n=1,4
               anabla2(i,n) = anabla2(i,n) + q(nb3,n)
  145       continue
         else
            nb3 = - nb3
            pp3 = ppg(nb3)
            do 150 n=1,4
               anabla2(i,n) = anabla2(i,n) + rghost(nb3,n)
  150       continue
         end if
c
c        Sum it all up.
c
c        anu(i) = dabs(pp1 + pp2 + pp3 - 3.0d0*pp(i))
c    &              / (pp1 + pp2 + pp3 + 3.0d0*pp(i))
c
         anu(i) = (dabs(pp1 - pp(i)) +
     &             dabs(pp2 - pp(i)) +
     &             dabs(pp3 - pp(i)))
     &             / (pp1 + pp2 + pp3 + 3.0d0*pp(i))
c
         do 160 n=1,4
            anabla2(i,n) = anabla2(i,n) - 3.0d0 * q(i,n)
  160    continue
  200 continue
c
c     Finally, loop over all real volumes and compute the artificial
c     dissipation operator.
c
      do 300 i=1,nvol
c
c        Loop over the faces of each control volume.
c
         do 290 k=1,3
            nb = neighbor(i,k)
            if(nb.lt.0) go to 250
c
c           If the neighbor is also a real volume, just do the
c           standard stuff.
c
            epsik2 = ak2 * max(anu(i),anu(nb))
            auxvar = ak4 - epsik2
            epsik4 = max(0.d0,auxvar)
c
            if(ivtstep.eq.1) then
               const4 = ((vol(i,1)/deltat(i)) + (vol(nb,1)/deltat(nb)))
     &                  / 2.0d0
            else
               const4 = (vol(i,1) + vol(nb,1)) / (2.0d0 * h)
            end if
            const2 = epsik2 * const4
            const4 = epsik4 * const4
            do 230 n=1,4
               dik2(n) = const2 * (q(nb,n) - q(i,n))
               dik4(n) = const4 * (anabla2(nb,n) - anabla2(i,n))
  230       continue
            go to 270
c
c           If the neighbor is a ghost volume, then apply the special
c           boundary procedure.
c
  250       continue
            nb = - nb
            epsik2 = ak2 * anu(i)
c
            if(ivtstep.eq.1) then
               const2 = epsik2 * (vol(i,1) + rghost(nb,5))
     &                  / (2.0d0 * deltat(i))
            else
               const2 = epsik2 * (vol(i,1) + rghost(nb,5)) / (2.0d0 * h)
            end if
            do 260 n=1,4
               dik2(n) = const2 * (rghost(nb,n) - q(i,n))
               dik4(n) = 0.0d0
  260       continue
c
c           Now, add these contributions to the artificial dissipation
c           operator.
c
  270       continue
            do 280 n=1,4
               dpen(i,n) = dpen(i,n) + (dik2(n) - dik4(n))
  280       continue
c
  290    continue
  300 continue
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine DISSMAVRI (Fred - 16-sep-2003)
c
c     This subroutine computes the dissipation operator for all the
c     "real" volumes using Mavriplis ideas and nomenclature.
c     (AIAA Journal, 1990)
c
c     Here I am taking full benefit of the edge based data, but there
c     are some IFs inside DOs due to boundary edges.
c     I'll fix that in the future...
c
      subroutine dissmavri
c
      include   'Common_rm_v2.0.f'
c
      integer   i,j, n1, n2, v1, v2
      dimension u(nmaxvl),v(nmaxvl),p(nmaxvl),
     &          ug(nmaxgh),vg(nmaxgh),pgh(nmaxgh),
     &          Ab(nmaxvl),Sp(nmaxvl),dw(nmaxvl,4),
     &          psens(nmaxvl),Spn(nmaxvl),Spd(nmaxvl)
      real*8    uav, vav, cav, delx, dely, expr, expr2,
     &          psmax, Abav, eps2b, eps1b, pdiff, psum
c
c     Zero out the dissipation array for upwind schemes
c
      if(ischeme.gt.2) then
         do j=1,4
            do i=1,nvol
               dpen(i,j) = 0.0d0
            enddo
         enddo
         return
      endif
c
c     Cells loop to initialize variables.
c
      do i=1,nvol
         u(i)    = q(i,2) / q(i,1)
         v(i)    = q(i,3) / q(i,1)
         p(i)    = (gamma-1.0d0) * (q(i,4)-0.5d0*q(i,1)*(u(i)*u(i)
     &                                                  +v(i)*v(i)))
         Ab(i)   = 0.0d0
         Sp(i)   = 0.0d0
         Spn(i)  = 0.0d0
         Spd(i)  = 0.0d0
         do j=1,4
            dw(i,j)   = -3.0d0 * q(i,j)
            dpen(i,j) =  0.0d0
         enddo
      enddo
      do i=1,nvolg
         ug(i)   = rghost(i,2) / rghost(i,1)
         vg(i)   = rghost(i,3) / rghost(i,1)
         pgh(i)   = (gamma-1.0d0) * (rghost(i,4)-0.5d0*rghost(i,1)*
     &                             (ug(i)*ug(i)+vg(i)*vg(i)))
      enddo
c
c     Edges loop to compute sums.
c
      do i=1,nsides
         n1 = iedges(i,1)
         n2 = iedges(i,2)
         v1 = iedges(i,3)
         v2 = iedges(i,4)
c
         if (v2.lt.0) then
            uav = 0.5d0 * (u(v1) + ug(-v2)) - vside(i,1)
            vav = 0.5d0 * (v(v1) + vg(-v2)) - vside(i,2)
            cav = 0.5d0 * (dsqrt(gamma*p(v1)/q(v1,1))+
     &                     dsqrt(gamma*pgh(-v2)/rghost(-v2,1)))
            delx = xynode(n2,1) - xynode(n1,1)
            dely = xynode(n2,2) - xynode(n1,2)
            expr = dabs(uav*dely-vav*delx)
     &             + cav*dsqrt(delx*delx+dely*dely)
            Ab(v1) = Ab(v1) + expr
c            Sp(v1) = Sp(v1) + pgh(-v2)
c-----modulo dentro do somatorio------------
            pdiff = dabs(pgh(-v2) - p(v1))
            psum  = pgh(-v2) + p(v1)
            Spn(v1) = Spn(v1) + pdiff
            Spd(v1) = Spd(v1) + psum
c-------------------------------------------
            do j=1,4
               dw(v1,j) = dw(v1,j) + rghost(-v2,j)
            enddo
         else
            uav = 0.5d0 * (u(v1) + u(v2)) - vside(i,1)
            vav = 0.5d0 * (v(v1) + v(v2)) - vside(i,2)
            cav = 0.5d0 * (dsqrt(gamma*p(v1)/q(v1,1))+
     &                     dsqrt(gamma*p(v2)/q(v2,1)))
            delx = xynode(n2,1) - xynode(n1,1)
            dely = xynode(n2,2) - xynode(n1,2)
            expr = dabs(uav*dely-vav*delx)
     &             + cav*dsqrt(delx*delx+dely*dely)
            Ab(v1) = Ab(v1) + expr
            Ab(v2) = Ab(v2) + expr
c            Sp(v1) = Sp(v1) + p(v2)
c            Sp(v2) = Sp(v2) + p(v1)
c-----modulo dentro do somatorio------------
            pdiff = dabs(p(v2) - p(v1))
            psum  = p(v2) + p(v1)
            Spn(v1) = Spn(v1) + pdiff
            Spd(v1) = Spd(v1) + psum
            Spn(v2) = Spn(v2) + pdiff
            Spd(v2) = Spn(v2) + psum
c-------------------------------------------
            do j=1,4
               dw(v1,j) = dw(v1,j) + q(v2,j)
               dw(v2,j) = dw(v2,j) + q(v1,j)
            enddo
         endif
c
      enddo
c
c     Cells loop to compute pressure sensor
c
      do i=1,nvol
c         psens(i) = dabs(Sp(i)-3.0d0*p(i))/(Sp(i)+3.0d0*p(i))
c-----modulo dentro do somatorio------------
         psens(i) = Spn(i)/Spd(i)
c-------------------------------------------
      enddo
c
c     Final edges loop
c
      do i=1,nsides
         v1 = iedges(i,3)
         v2 = iedges(i,4)
         if (v2.lt.0) then
            Abav = Ab(v1)
c-----ponderacao do jameson-----------------
c            Abav = vol(v1,1)/deltat(v1)
c            Abav = vol(v1,1)/h
c-------------------------------------------
            psmax = psens(v1)
c-----sensor de pressao do Stolcis------------------------
c            psmax = dabs(pgh(-v2)-p(v1))/(pgh(-v2)+p(v1))
c---------------------------------------------------------
            eps1b = ak2 * psmax
            do j=1,4
               expr2 = Abav * eps1b * (rghost(-v2,j)-q(v1,j))
               dpen(v1,j) = dpen(v1,j) + expr2
            enddo
         else
            Abav = 0.5d0 * (Ab(v1)+Ab(v2))
c-----ponderacao do jameson------------------------------------------------
c            Abav = 0.5d0 * (vol(v1,1)/deltat(v1) + vol(v2,1)/deltat(v2))
c            Abav = 0.5d0 * (vol(v1,1)/h + vol(v2,1)/h)
c--------------------------------------------------------------------------
            psmax = dmax1(psens(v1),psens(v2))
c-----sensor de pressao do Stolcis------------------------
c            psmax = dabs(p(v2)-p(v1))/(p(v2)+p(v1))
c---------------------------------------------------------
            eps1b = ak2 * psmax
            eps2b = dmax1(0.0d0,(ak4-eps1b))
            do j=1,4
               expr2 = Abav * (eps1b * (q(v2,j)-q(v1,j))
     &                       - eps2b * (dw(v2,j)-dw(v1,j)))
               dpen(v1,j) = dpen(v1,j) + expr2
               dpen(v2,j) = dpen(v2,j) - expr2
            enddo
         endif
      enddo
c
      return
      end
c ---------------------------------------------------------------------
c
c     Subroutine FLUX
c
c     This subroutine computes the E and F flux vectors at cell
c     interfaces, using the averaged value of the conserved properties
c     at that interface. The computation is performed for one single
c     interface each time the routine is called.
c
c     This routine may have to be modified in the future such that it
c     could compute all forms of flux vectors that may be required in
c     the upwind schemes.
c
      subroutine flux
c
      include 'Common_rm_v2.0.f'
c
c
c     Perform preliminary computations.
c
      const = gamma - 1.0d0
      u = qav(2) / qav(1)
      v = qav(3) / qav(1)
      p = const * (qav(4) - 0.5d0*(qav(2)*qav(2)
     &                           + qav(3)*qav(3)) / qav(1))
c
c     Compute the contravariant velocities
c
      uu = u - uedge
      vv = v - vedge
c
c     Compute the E and F flux vectors.
c
      evec(1) = qav(1)*uu
      evec(2) = qav(2)*uu + p
      evec(3) = qav(3)*uu
      evec(4) = (qav(4) + p) * uu + uedge*p
      fvec(1) = qav(1)*vv
      fvec(2) = qav(2)*vv
      fvec(3) = qav(3)*vv + p
      fvec(4) = (qav(4) + p) * vv + vedge*p
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine INDAT
c
c     This subroutine reads all the input data, writes an echo of this
c     input data for later verification, and calls the appropriate
c     routines for cell "volume" computation and for the set up of
c     initial conditions for the present run.
c
      subroutine indat
c
      include 'Common_rm_v2.0.f'
c
c     Temporary variable declaration
c
      real*4 pretime,timecomp(2)
c
c     Read the input data (unit = 1).
c


       
      read(1,901) nstart,nmiter,nsave,nhist
      read(1,902) nnodes,nvol
      read(1,903) nvisc
      read(1,906) 
      read(1,903) iflow
      read(1,907) 
      read(1,903) iflow2
      read(1,906) 
      read(1,903) ivtstep
      read(1,906) 
      read(1,903) itime
      read(1,906) 
      read(1,904) h
      read(1,903) ischeme
      read(1,909) 
      read(1,903) idissip
      read(1,906) 
      read(1,905) ak2,ak4
      read(1,903) irsmth
      read(1,906) 
      read(1,910) ljacrs,epsirs
      read(1,903) irstag
      read(1,906) 
      read(1,903) ibcchr
      read(1,906) 
      read(1,904) pg
      read(1,903) nitj
      read(1,903) ivmesh
      read(1,906) 
      read(1,903) igcl
      read(1,906) 
      read(1,911) fsmach,rey,gamma,pexit,xea
      read(1,903) imotion
      read(1,907)
      read(1,903) idof
      read(1,906) 
      read(1,903) imove
      read(1,908)     
      read(1,911) alpha,akc,alpham,ymax,tmax      
      read(1,912) epscon,epsblw
      read(1,903) iresp
      read(1,907)      
      read(1,914) nintsave
      read(1,*)
      read(1,903) ipress
      read(1,907)
      read(1,914) npressave
      close(1)
c
c     Write this data for posterity.
c
      write(10,921) nstart,nmiter,nsave,nhist
      write(10,999)
      write(10,922) nnodes,nvol
      write(10,999)
      write(10,923) nvisc
      write(10,999)
      write(10,924) iflow
      write(10,999)
      write(10,925) iflow2
      write(10,999)
      write(10,926) ivtstep
      write(10,999)
      write(10,927) itime
      write(10,999)
      write(10,928) h
      write(10,999)
      write(10,929) ischeme
      write(10,999)
      write(10,930) idissip
      write(10,999)
      write(10,931) ak2,ak4
      write(10,999)
      write(10,932) irsmth
      write(10,999)
      write(10,933) ljacrs,epsirs
      write(10,999)
      write(10,934) irstag
      write(10,999)
      write(10,935) ibcchr
      write(10,999)
      write(10,936) pg,nitj
      write(10,999)
      write(10,937) ivmesh
      write(10,999)
      write(10,938) igcl
      write(10,999)
      write(10,939) fsmach,rey,gamma,pexit,xea
      write(10,999)
      write(10,948) imotion
      write(10,999)
      write(10,940) idof
      write(10,999)
      write(10,941) imove
      write(10,999)
      write(10,942) alpha,akc,alpham,ymax,tmax
      write(10,999)
      write(10,943) epscon,epsblw
      write(10,999)
      write(10,944) iresp
      write(10,999)
      write(10,945) nintsave
      write(10,999)
      write(10,946) ipress
      write(10,999)
      write(10,947) npressave
      close(10)
c
c     Make things compatible with the code needs.
c
      pi = dacos(-1.0d0)
c
      alpha  = alpha * pi / 180.0d0
      alpham = alpham * pi / 180.0d0
c     
c     Read the mesh and compute the "volume" of all triangular control
c     volumes.
c
      call meshio
c
c     Here comes the first time it writes the elapsed time.
c  
        print*, "saiu mesh io"
      pretime = dtime(timecomp)
      if(nstart.eq.0) then
         open(22,file='Fort.22')
         write(22,920) pretime
         close(22)
      endif
c
c     Start the edge-based database. Eventually, this task will be
c     performed by the grid generation routine in the future.
c
      call preproc
c
        print*, "saiu preproc"
c     Initialize the flow variables.
c
      call start
c
c     If we are running a 2nd order upwind scheme option, then one has
c     to identify the triangles which will be used for reconstruction
c     for each edge. As with the call above, eventually this may be
c     done by the grid generation program.
c
      if(ischeme.ge.5) call preproc2
c
      return
c
c     Formats.
c
  901 format(3(36x,i10,/),/,36x,i10,/)
  902 format(2(36x,i10,/))
  903 format(36x,i1,/)
  904 format(36x,f10.8,/)
  905 format(2(36x,f10.8,/))
  906 format(/)
  907 format(//)
  908 format(///)
  909 format(////)
  910 format(36x,i2,/,36x,f10.8,/)
  911 format(5(36x,f10.8,/))
  912 format(2(36x,e11.5,/))
  913 format(4(36x,f10.8,/),36x,f10.5,/)
  914 format(36x,i10)
  920 format(1x,'Preprocessing time (CPU s):',1x,e9.3)
  921 format(//,15x,'INPUT DATA',///,
     &       'initial iteration number          = ',i10,/,
     &       'additional iteration number       = ',i10,/,
     &       'saving frequency for restart      = ',i10,/,
     &       'saving frequency for',/,
     &       'aerodynamic coefficients history  = ',i10)
  999 format('-------------------------------------------------')
  922 format('number of nodes in the mesh       = ',i10,/,
     &       'number of volumes in the mesh     = ',i10)
  923 format('physical model                    = ',i1,/,
     &       '(0) Euler',/,
     &       '(1) Navier-Stokes')
  924 format('flow situation                    = ',i1,/,
     &       '(1) external',/,
     &       '(2) internal, supersonic entrance',/,
     &       '(3) internal, subsonic entrance')
  925 format('flow type                         = ',i1,/,
     &       '(0) steady',/,
     &       '(1) unsteady')
  926 format('time-marching method              = ',i1,/,
     &       '(0) fixed',/,
     &       '(1) variable')
  927 format('time-marching scheme              = ',i1,/,
     &       '(0) 5-step Runge-Kutta',/,
     &       '(1) 4-step Runge-Kutta')
  928 format('time step (or CFL number)         = ',f10.8)
  929 format('space discretization              = ',i1,/,
     &       '(1) Jameson',/,
     &       '(2) MacCormack',/,
     &       '(3) 1st-order Van Leer FVS',/,
     &       '(4) Liou-Steffen FVS',/,
     &       '(5) 2nd-order Van Leer FVS')
  930 format('artificial dissipation            = ',i1,/,
     &       '(0) Jameson',/,
     &       '(1) Mavriplis')
  931 format('ak2 constant                      = ',f10.8,/,
     &       'ak4 constant                      = ',f10.8)
  932 format('implicit residual smoothing       = ',i1,/,
     &       '(0) off',/,
     &       '(1) on')
  933 format('number of Jacobi iterations       = ',i2,/,
     &       'irs epsilon                       = ',f10.8)
  934 format('irs effectiveness                 = ',i1,/,
     &       '(0) alternate step',/,
     &       '(1) every step')
  935 format('Entrance/exit boundary conditions = ',i1,/,
     &       '(0) zero-order',/,
     &       '(1) characteristic relations')
  936 format('grid stifness factor              = ',f10.8,/,
     &       'number of Jacobi',/,
     &       'iterations to move the mesh       = ',i1)
  937 format('node veloctiy calculation         = ',i1,/,
     &       '(1) 1st-order',/,
     &       '(2) 2nd-order')
  938 format('area calculation                  = ',i1,/,
     &       '(0) conventional',/,
     &       '(1) geometric conservation law')
  939 format('Mach number                       = ',f10.8,/,
     &       'Reynolds number                   = ',f10.8,/,
     &       'ratio of specific heats           = ',f10.8,/,
     &       'exit pressure                     = ',f10.8,/,
     &       'position of elastic axis          = ',f10.8)
  948 format('motion type                       = ',i1,/,
     &       '(0) altogether',/,
     &       '(1) displacement only',/,
     &       '(2) velocity only')
  940 format('mode                              = ',i1,/,
     &       '(0) plunging',/,
     &       '(1) pitching')
  941 format('motion kind                       = ',i1,/,
     &       '(0) harmonic oscillation',/,
     &       '(1) exponentially-shaped pulse',/,
     &       '(2) discrete step',/,
     &       '(3) unit sample')
  942 format('(initial) angle of attack         = ',f10.8,/,
     &       'reduced frequency                 = ',f10.8,/,
     &       'maximum angle of attack           = ',f10.8,/,
     &       'maximum translation               = ',f10.8,/,
     &       'pulse duration                    = ',f10.5)
  943 format('convergence residue               = ',e11.5,/,
     &       'divergence residue                = ',e11.5)
  944 format('save intermediary results         = ',i1,/,
     &       '(0) no save',/,
     &       '(1) save')
  945 format('saving frequency for',/,
     &       'intermediary results              = ',i10)
  946 format('save pressure distributions       = ',i1,/,
     &       '(0) no save',/,
     &       '(1) save')
  947 format('saving frequency for',/,
     &       'pressure distributions            = ',i10)
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine MAXR
c
c     This subroutine computes the maximum residue in the field.
c     We will be testing here the (delta_rho / rho) as the
c     representative residue quantity. Then, this is like a
c     L_infty norm of the residue (it is not a L_2 norm).
c
      subroutine maxr
c
      include 'Common_rm_v2.0.f'
c
c     Compute the representative quantity for the whole field.
c
      do 20 i=1,nvol
         qaux(i,1) = dabs(rhs(i,1) / q(i,1))
   20 continue
c
c     Initialize variables for search.
c
      irhsmax = 0
      rhsmax = 0.0d0
c
c     Perform search of maximum representative quantity in the field.
c
      do 60 i=1,nvol
         if(qaux(i,1).gt.rhsmax) then
            rhsmax = qaux(i,1)
            irhsmax = i
         end if
   60 continue
c
c     Compute log_10 of maximum residue and print information.
c
      aux1 = dlog10(rhsmax)
c
      pi = acos(-1.d0)
      alphai = alphat*180.d0/pi
c
      write(*,900) niter,alphai,aux1,rhsmax,irhsmax
      write(9,901) niter,aux1,rhsmax,irhsmax
c
c     Compute, now, the L_2 norm of the residue. This will still use
c     the representative quantity previously defined, i.e., the density
c     residue normalized by the density itself.
c
      sum = 0.0d0
      do 110 i=1,nvol
         sum = sum + qaux(i,1)*qaux(i,1)
  110 continue
c
      rdl2n = dsqrt(sum/(dble(nvol)))
      aux1 = dlog10(rdl2n)
      write(11,911) niter,aux1,rdl2n
c
      return
c
c     Formats.
c
  900 format(1x,i6,2x,f6.3,2x,f13.7,5x,e15.7,5x,i6)
  901 format(1x,i6,2x,f13.7,5x,e15.7,5x,i6)
  911 format(1x,i6,2x,f13.7,5x,e15.7)
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine MESHIO
c
c     This subroutine will read the grid file and the connectivity
c     tables. It will also compute the "volume" of all triangles.
c
      subroutine meshio
c
      include 'Common_rm_v2.0.f'
c
c     Read grid file (unit = 2).
c
      read(2,911) itest1
      do 20 i=1,nnodes
         print*, "2 - point ", i
         read(2,905) ii,xynode(i,1),xynode(i,2)
c
         xhist(i,1) = xynode(i,1)
         yhist(i,1) = xynode(i,2)
         x0(i) = xynode(i,1)
         y0(i) = xynode(i,2)
c
   20 continue
      close(2)
c
c     Read connectivity table (unit = 3).
c
      read(3,911) itest1
      if(itest1.ne.nvol) then
         write(*,912) itest1
         stop
      end if
c
      do 40 i=1,nvol
         print*, "3 - panel ", i
         read(3,915) ii,(itable(i,n),n=1,3)
   40 continue
      close(3)
c
c     Read table with the neighbors of each control volume (unit = 4).
c
      do 60 i=1,nvol
         print*, "4 - neighbours ", i
         read(4,921) ii,(neighbor(i,n),n=1,3)
   60 continue
      close(4)
c
c     Read the ghost volume information data (unit = 5).
c
      read(5,931) nvolg
      do 80 i=1,nvolg
         print*, "5 - ghosts ", i
         read(5,932) ii,(ighost(i,n),n=1,4)
   80 continue
      close(5)
c
c     For the unsteady flow case, we must generate the connectivity
c     table of the nodes.
c
      if (iflow2.eq.1) then
         call gen_f13
      endif
c
c     Compute "volume" (or the area, if you prefer) of each triangle.
c
        print*, "calculando volumes"
      do 140 i=1,nvol
         ina = itable(i,1)
         inb = itable(i,2)
         inc = itable(i,3)
c
         ax = xynode(ina,1)
         ay = xynode(ina,2)
         bx = xynode(inb,1)
         by = xynode(inb,2)
         cx = xynode(inc,1)
         cy = xynode(inc,2)
c
         call volumes(ax,ay,bx,by,cx,cy,area)
         vol(i,1) = area
         vol(i,2) = vol(i,1)
  140 continue
c
c     Now, fill the appropriate portion of the "rghost" array with
c     the "volume" (or "area") of that ghost volume.
c
        print*, "volumes calculados"
      do 160 ig=1,nvolg
         intngb = ighost(ig,2)
         rghost(ig,5) = vol(intngb,1)
  160 continue
c
c     For variable time-stepping, compute the smallest characteristic
c     length of each real (or interior) control volume. The idea here
c     is that the smallest side of the triangle, or the smallest
c     distance between the centroid of the present triangle and that
c     of its neighbors (whichever is smaller), will be taken as the
c     characteristic length used for the space-varying time-step
c     calculation for that cell (based on a constant CFL number).
c
      if(ivtstep.eq.1) then
         do 250 i=1,nvol
c
c           Find smallest side of the current triangle.
c
            in1 = itable(i,1)
            in2 = itable(i,2)
            in3 = itable(i,3)
            it1 = in1
            it2 = in2
c
            auxx = xynode(it2,1) - xynode(it1,1)
            auxy = xynode(it2,2) - xynode(it1,2)
            sidemin = dsqrt(auxx*auxx + auxy*auxy)
c
            it1 = in2
            it2 = in3
            do 220 k=2,3
               auxx = xynode(it2,1) - xynode(it1,1)
               auxy = xynode(it2,2) - xynode(it1,2)
               sidecur = dsqrt(auxx*auxx + auxy*auxy)
               if(sidecur.lt.sidemin) sidemin = sidecur
               if(k.eq.2) then
                  it1 = in3
                  it2 = in1
               end if
  220       continue
c
c           Find smallest distance between the centroid of the present
c           triangle and that of its neighbors. If the neighbor is a
c           ghost volume, then it is not considered in the present
c           calculation.
c
            distmin = 100.0d0 * sidemin
            centrox = (xynode(in1,1) + xynode(in2,1) + xynode(in3,1))
     &                / 3.0d0
            centroy = (xynode(in1,2) + xynode(in2,2) + xynode(in3,2))
     &                / 3.0d0
c
            do 240 k=1,3
               nb = neighbor(i,k)
               if(nb.lt.0) go to 230
               it1 = itable(nb,1)
               it2 = itable(nb,2)
               it3 = itable(nb,3)
               centnbx = (xynode(it1,1) + xynode(it2,1) +
     &                    xynode(it3,1)) / 3.0d0
               centnby = (xynode(it1,2) + xynode(it2,2) +
     &                    xynode(it3,2)) / 3.0d0
               auxx = centnbx - centrox
               auxy = centnby - centroy
               dist = dsqrt(auxx*auxx + auxy*auxy)
               if(dist.lt.distmin) distmin = dist
  230          continue
  240       continue
c
c           Finally, pick the smallest between the two results as the
c           desired characteristic length.
c
            if(sidemin.lt.distmin) then
               clength(i) = sidemin
            else
               clength(i) = distmin
            end if
  250    continue
      end if
c
      return
c
c     Formats.
c
  905 format(i6,2x,2(e16.8,2x))
  911 format(i7)
  912 format(///,1x,'MESHIO, number of triangles in the connectivity',
     & ' table (',i7,')',/,1x,'does not agree with the original input')
  915 format((i7,2x),3(i6,2x))
  921 format(4(i7,2x))
  931 format(i7)
  932 format(5(i7,2x))
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine PREPROC
c
c     This subroutine creates the edge-based data structure using the
c     data that the code is usually reading.
c     This essentially amounts to create a new form of indexing, or
c     connectivity, which will be used for the flux calculations (I
c     mean, both convective and dissipative fluxes and the residue
c     calculations).
c     At this point in time, there will be no attempt to complete
c     suppress the other forms of indexing available in the code.
c     This means that there are going to be more arrays (and more
c     memory usage) than it is strictly necessary. But, this is
c     something to fix in the future...
c
      subroutine preproc
c
      include 'Common_rm_v2.0.f'
c
c     I will not try to be fancy here. Just loop over the control
c     volumes and create the new table.
c
      nsides = 0
      do 400 i=1,nvol
      print*, "preproc ", i
c
c        For each control volume, loop over each side.
c
         do 350 k=1,3
c
c           For each side, get the necessary connectivity information.
c
            if(k.eq.1) then
               in1 = itable(i,1)
               in2 = itable(i,2)
               nb = neighbor(i,1)
            else
               if(k.eq.2) then
                  in1 = itable(i,2)
                  in2 = itable(i,3)
                  nb = neighbor(i,2)
               else
                  in1 = itable(i,3)
                  in2 = itable(i,1)
                  nb = neighbor(i,3)
               end if
            end if
c
c           Now, check whether the edge already exists.
c
            if(nsides.gt.0) then
               do 100 j=1,nsides
                  itry1 = iedges(j,1)
                  itry2 = iedges(j,2)
                  if((in1.eq.itry2).and.(in2.eq.itry1).and.(nb.gt.0))
     &              go to 340
  100          continue
            end if
c
c           If we get to this point, the edge does not exist yet.
c           Hence, add it up to the list of edges.
c
            nsides = nsides + 1
            iedges(nsides,1) = in1
            iedges(nsides,2) = in2
            iedges(nsides,3) = i
            iedges(nsides,4) = nb
c
  340       continue
  350    continue
  400 continue
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine RESIDUE
c
c     This subroutine computes the flux vectors for all interior
c     volumes, the flux vectors for the ghost volumes, the convective
c     and the dissipative operators to form the RHS, and finally it
c     determines the maximum residue in the whole field for convergence
c     test purposes.
c
      subroutine residue
c
      include 'Common_rm_v2.0.f'
c
c     Compute convective operator.
c
      call convec
c
c     If the variable time step option is selected, we must compute
c     the delta_t's here, because Dissip will need them.
c
      if(ivtstep.eq.1) call vartime
c
c     Compute dissipative operator.
c
      if (idissip.eq.0) then
         call dissip
      elseif (idissip.eq.1) then
         call dissmavri
      endif
c
c     Compute current residue.
c
      if (iflow2.eq.1 .and. niter.gt.1) then
         do 100 i=1,nvol
            rhs(i,1) = (cpen(i,1) - dpen(i,1)) / vol(i,2)
            rhs(i,2) = (cpen(i,2) - dpen(i,2)) / vol(i,2)
            rhs(i,3) = (cpen(i,3) - dpen(i,3)) / vol(i,2)
            rhs(i,4) = (cpen(i,4) - dpen(i,4)) / vol(i,2)
  100    continue
      else
        do 110 i=1,nvol
           rhs(i,1) = (cpen(i,1) - dpen(i,1)) / vol(i,1)
           rhs(i,2) = (cpen(i,2) - dpen(i,2)) / vol(i,1)
           rhs(i,3) = (cpen(i,3) - dpen(i,3)) / vol(i,1)
           rhs(i,4) = (cpen(i,4) - dpen(i,4)) / vol(i,1)
  110   continue
      endif
c
c     Compute maximum residue in the field (unFortunately, before
c     that, I have to eliminate the unpleasant case of the first
c     iteration when iflow=3, i.e., internal flow with subsonic
c     entrance).
c
      if(iflow.eq.3.and.niter.eq.1) then
         rhsmax = 1.0d0
         return
      end if
c
      call maxr
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine RESMTH
c
c     This subroutine performs the implicit residual smoothing (or
c     implicit residual averaging). The current residue, called R(Q)
c     in my notes, is the input data, and the smoothed residue,
c     called RBAR(Q) in my notes, is the output. The solution is
c     achieved by Jacobi iterations, and the number of iterations
c     performed is specified by the parameter "ljacrs".
c
      subroutine resmth
c
      include 'Common_rm_v2.0.f'
c
c
      dimension rbar(nmaxvl,4),rblm1(nmaxvl,4)
c
c     Some initial definitions.
c
      const = 1.0d0 / (1.0d0 + 3.0d0*epsirs)
      lastm1 = ljacrs - 1
c
c     Create initial guess for the smoothed residual and zero out
c     the "rbar" array.
c
      do 50 i=1,nvol
         rblm1(i,1) = rhs(i,1)
         rblm1(i,2) = rhs(i,2)
         rblm1(i,3) = rhs(i,3)
         rblm1(i,4) = rhs(i,4)
         rbar(i,1) = 0.0d0
         rbar(i,2) = 0.0d0
         rbar(i,3) = 0.0d0
         rbar(i,4) = 0.0d0
   50 continue
c
c     Perform (ljacrs - 1) Jacobi iterations and updates.
c
      do 150 l=1,lastm1
c
c        Initially, loop over all volumes and add up the contribution
c        of their neighbors. (Observe that this loop probably will
c        not vectorize.)
c
         do 100 i=1,nvol
            do 80 k=1,3
               nb = neighbor(i,k)
               if(nb.lt.0) go to 70
               rbar(i,1) = rbar(i,1) + rblm1(nb,1)
               rbar(i,2) = rbar(i,2) + rblm1(nb,2)
               rbar(i,3) = rbar(i,3) + rblm1(nb,3)
               rbar(i,4) = rbar(i,4) + rblm1(nb,4)
   70          continue
   80       continue
  100    continue
c
c        Compute new residues for all cells.
c
         do 120 i=1,nvol
            rbar(i,1) = const * (rhs(i,1) + epsirs*rbar(i,1))
            rbar(i,2) = const * (rhs(i,2) + epsirs*rbar(i,2))
            rbar(i,3) = const * (rhs(i,3) + epsirs*rbar(i,3))
            rbar(i,4) = const * (rhs(i,4) + epsirs*rbar(i,4))
  120    continue
c
c        Now, perform the update.
c
         do 130 i=1,nvol
            rblm1(i,1) = rbar(i,1)
            rblm1(i,2) = rbar(i,2)
            rblm1(i,3) = rbar(i,3)
            rblm1(i,4) = rbar(i,4)
  130    continue
c
c        And, zero the "rbar" array for the next iteration.
c
         do 140 i=1,nvol
            rbar(i,1) = 0.0d0
            rbar(i,2) = 0.0d0
            rbar(i,3) = 0.0d0
            rbar(i,4) = 0.0d0
  140    continue
  150 continue
c
c     Perform the last Jacobian iteration.
c
c     As before, initially loop over all the cells and accumulate
c     the contribution of the neighbors.
c
      do 300 i=1,nvol
         do 280 k=1,3
            nb = neighbor(i,k)
            if(nb.lt.0) go to 270
            rbar(i,1) = rbar(i,1) + rblm1(nb,1)
            rbar(i,2) = rbar(i,2) + rblm1(nb,2)
            rbar(i,3) = rbar(i,3) + rblm1(nb,3)
            rbar(i,4) = rbar(i,4) + rblm1(nb,4)
  270       continue
  280    continue
  300 continue
c
c     Then, compute the new residue for all cell.
c
      do 350 i=1,nvol
         rbar(i,1) = const * (rhs(i,1) + epsirs*rbar(i,1))
         rbar(i,2) = const * (rhs(i,2) + epsirs*rbar(i,2))
         rbar(i,3) = const * (rhs(i,3) + epsirs*rbar(i,3))
         rbar(i,4) = const * (rhs(i,4) + epsirs*rbar(i,4))
  350 continue
c
c     Finally, update the true RHS array and return.
c
      do 400 i=1,nvol
         rhs(i,1) = rbar(i,1)
         rhs(i,2) = rbar(i,2)
         rhs(i,3) = rbar(i,3)
         rhs(i,4) = rbar(i,4)
  400 continue
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine RK4S
c
c     This subroutine performs one time step according to the 4-stage,
c     2nd-order, Runge-Kutta time-stepping scheme described by
c     Jameson, Mavriplis and co-workers (this is not meant to be a
c     true reference list). The alpha_i parameters are set directly
c     in the beginning of the subroutine. So, there was no attempt
c     to make this a "general" Runge-Kutta time-stepping routine.
c
      subroutine rk4s
c
      include 'Common_rm_v2.0.f'
c
c
      dimension alphp(4),const1(nmaxvl),const2(nmaxvl)
c
c     Definition of alpha_i parameters.
c
      alphp(1) = 1.0d0 / 4.0d0
      alphp(2) = 1.0d0 / 3.0d0
      alphp(3) = 1.0d0 / 2.0d0
      alphp(4) = 1.0d0
c
c     First, we have to save the current solution in the auxiliary
c     array "qaux" (i.e., qaux = q(0)).
c
      do 10 i=1,nvol
         qaux(i,1) = q(i,1)
         qaux(i,2) = q(i,2)
         qaux(i,3) = q(i,3)
         qaux(i,4) = q(i,4)
   10 continue
c
c     If the variable time step option is in effect, things should be
c     done in a different form.
c
      if(ivtstep.eq.1) go to 500
c
c     OK, if you get to this point, then this is a constant (or global)
c     time step case.
c
c     The next decision which must be made is whether there is implicit
c     residual smoothing or not.
c
      if(irsmth.eq.1) go to 300
c
c     This is the case without residual smoothing, then.
c
c     The quantity that comes into this routine is already the current
c     residue. Therefore, we can go straight to stage # 1.
c
c     Stage # 1.
c
      const = h * alphp(1)
c                                
c     Here we are considering the case where we have unsteady flow.
c
      if (iflow2.eq.1 .and. niter.gt.1) then
        do 20 i=1,nvol
          const2(i) = vol(i,1)/vol(i,2) 
   20   continue
      else
        do 30 i=1,nvol              
          const2(i) = 1.d0
   30   continue
      endif
c
      do 40 i=1,nvol
         q(i,1) = const2(i)*qaux(i,1) - const*rhs(i,1)
         q(i,2) = const2(i)*qaux(i,2) - const*rhs(i,2)
         q(i,3) = const2(i)*qaux(i,3) - const*rhs(i,3)
         q(i,4) = const2(i)*qaux(i,4) - const*rhs(i,4)
   40 continue
c
c     Now, do Stages # 2 to 4.
c
      do 200 k=2,4
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
         call convec
c
c        Perform the Runge-Kutta stage.
c
         const = h * alphp(k)
c
         if (iflow2.eq.1 .and. niter.gt.1) then
           do 50 i=1,nvol
             const1(i) = const/vol(i,2)                  
   50      continue
         else
           do 60 i=1,nvol              
             const1(i) = const/vol(i,1)
   60      continue
         endif
c
         do 150 i=1,nvol
            c1 = const2(i)*qaux(i,1)
            c2 = const2(i)*qaux(i,2)
            c3 = const2(i)*qaux(i,3)
            c4 = const2(i)*qaux(i,4)
c 
            q(i,1) = c1 - const1(i)*(cpen(i,1)-dpen(i,1))
            q(i,2) = c2 - const1(i)*(cpen(i,2)-dpen(i,2))
            q(i,3) = c3 - const1(i)*(cpen(i,3)-dpen(i,3))
            q(i,4) = c4 - const1(i)*(cpen(i,4)-dpen(i,4))
  150    continue
  200 continue
c
      return
c
c     This is the case with implicit residual smoothing (or implicit
c     residual averaging) but without variable time step.
c
  300 continue
c
c     The quantity that comes into this routine is the current
c     unsmoothed residue. Hence, we must smooth it first.
c
      call resmth
c
c     Stage # 1.
c
      const = h * alphp(1)
      do 350 i=1,nvol
         q(i,1) = qaux(i,1) - const * rhs(i,1)
         q(i,2) = qaux(i,2) - const * rhs(i,2)
         q(i,3) = qaux(i,3) - const * rhs(i,3)
         q(i,4) = qaux(i,4) - const * rhs(i,4)
  350 continue
c
c     Now, do Stages # 2 to 5.
c
      do 400 k=2,4
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
         call convec
c
c        Compute the new residue (unsmoothed).
c
         do 370 i=1,nvol
            rhs(i,1) = (cpen(i,1)-dpen(i,1)) / vol(i,1)
            rhs(i,2) = (cpen(i,2)-dpen(i,2)) / vol(i,1)
            rhs(i,3) = (cpen(i,3)-dpen(i,3)) / vol(i,1)
            rhs(i,4) = (cpen(i,4)-dpen(i,4)) / vol(i,1)
  370    continue
c
c        Smooth it when appropriate.
c
         if(k.eq.4) then
            call resmth
         else
            if(irstag.eq.1) call resmth
         end if
c
c        Perform the Runge-Kutta stage.
c
         const = h * alphp(k)
         do 390 i=1,nvol
            q(i,1) = qaux(i,1) - const * rhs(i,1)
            q(i,2) = qaux(i,2) - const * rhs(i,2)
            q(i,3) = qaux(i,3) - const * rhs(i,3)
            q(i,4) = qaux(i,4) - const * rhs(i,4)
  390    continue
  400 continue
c
      return
c
c     This is the variable time step case.
c
c     If this option is selected, remember that "h" is the CFL number.
c
  500 continue
c
c     Here, we must also decide whether there is going to be residual
c     smoothing or not.
c
      if(irsmth.eq.1) go to 800
c
c     OK, then it is going to be variable time stepping but without
c     residual smoothing.
c
c     The quantity that comes into this routine is already the current
c     residue. Therefore, we can go straight to stage # 1.
c
c     Stage # 1.
c
      do 550 i=1,nvol
         const = alphp(1) * deltat(i)
         q(i,1) = qaux(i,1) - const*rhs(i,1)
         q(i,2) = qaux(i,2) - const*rhs(i,2)
         q(i,3) = qaux(i,3) - const*rhs(i,3)
         q(i,4) = qaux(i,4) - const*rhs(i,4)
  550 continue
c
c     Now, do Stages # 2 to 5.
c
      do 700 k=2,4
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
         call convec
c
c        Perform the Runge-Kutta stage.
c
         do 650 i=1,nvol
            const = alphp(k) * deltat(i) / vol(i,1)
            q(i,1) = qaux(i,1) - const*(cpen(i,1)-dpen(i,1))
            q(i,2) = qaux(i,2) - const*(cpen(i,2)-dpen(i,2))
            q(i,3) = qaux(i,3) - const*(cpen(i,3)-dpen(i,3))
            q(i,4) = qaux(i,4) - const*(cpen(i,4)-dpen(i,4))
  650    continue
  700 continue
c
      return
c
c     Finally, this is the case with variable time stepping and with
c     implicit residual smoothing (or implicit residual averaging).
c
  800 continue
c
c     The quantity that comes into this routine is the current
c     unsmoothed residue. Hence, the first step is to smooth it,
c     before performing the first stage of the R-K time stepping.
c
      call resmth
c
c     Stage # 1.
c
      do 850 i=1,nvol
         const = alphp(1) * deltat(i)
         q(i,1) = qaux(i,1) - const*rhs(i,1)
         q(i,2) = qaux(i,2) - const*rhs(i,2)
         q(i,3) = qaux(i,3) - const*rhs(i,3)
         q(i,4) = qaux(i,4) - const*rhs(i,4)
  850 continue
c
c     Now, perform Stages # 2 to 5.
c
      do 1000 k=2,4
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
         call convec
c
c        Compute the new (unsmoothed) residue using the current fluxes.
c
         do 900 i=1,nvol
            rhs(i,1) = (cpen(i,1) - dpen(i,1)) / vol(i,1)
            rhs(i,2) = (cpen(i,2) - dpen(i,2)) / vol(i,1)
            rhs(i,3) = (cpen(i,3) - dpen(i,3)) / vol(i,1)
            rhs(i,4) = (cpen(i,4) - dpen(i,4)) / vol(i,1)
  900    continue
c
c        Perform the appropriate residual smoothing for this stage.
c
         if(k.eq.4) then
            call resmth
         else
            if(irstag.eq.1) call resmth
         end if
c
c        Perform the Runge-Kutta stage.
c
         do 950 i=1,nvol
            const = alphp(k) * deltat(i)
            q(i,1) = qaux(i,1) - const * rhs(i,1)
            q(i,2) = qaux(i,2) - const * rhs(i,2)
            q(i,3) = qaux(i,3) - const * rhs(i,3)
            q(i,4) = qaux(i,4) - const * rhs(i,4)
  950    continue
 1000 continue
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine RK5S
c
c     This subroutine performs one time step according to the 5-stage,
c     2nd-order, Runge-Kutta time-stepping scheme described by
c     Jameson, Mavriplis and co-workers (this is not meant to be a
c     true reference list). The alpha_i parameters are set directly
c     in the beginning of the subroutine. So, there was no attempt
c     to make this a "general" Runge-Kutta time-stepping routine.
c
      subroutine rk5s
c
      include 'Common_rm_v2.0.f'
c
      dimension alphp(5),const1(nmaxvl),const2(nmaxvl)
c
c     Definition of alpha_i parameters.
c
      alphp(1) = 1.0d0 / 4.0d0
      alphp(2) = 1.0d0 / 6.0d0
      alphp(3) = 3.0d0 / 8.0d0
      alphp(4) = 1.0d0 / 2.0d0
      alphp(5) = 1.0d0
c
c     First, we have to save the current solution in the auxiliary
c     array "qaux" (i.e., qaux = q(0)).
c
      do 10 i=1,nvol
         qaux(i,1) = q(i,1)
         qaux(i,2) = q(i,2)
         qaux(i,3) = q(i,3)
         qaux(i,4) = q(i,4)
   10 continue
c
c     If the variable time step option is in effect, things should be
c     done in a different form.
c
      if(ivtstep.eq.1) go to 500
c
c     OK, if you get to this point, then this is a constant (or global)
c     time step case.
c
c     The next decision which must be made is whether there is implicit
c     residual smoothing or not.
c
      if(irsmth.eq.1) go to 300
c
c     This is the case without residual smoothing, then.
c
c     The quantity that comes into this routine is already the current
c     residue. Therefore, we can go straight to stage # 1.
c
c     Stage # 1.
c
      const = h * alphp(1)
c                                
c     Here we are considering the case were we have unsteady flow.
c
      if (iflow2.eq.1 .and. niter.gt.1) then
        do 20 i=1,nvol
          const2(i) = vol(i,1)/vol(i,2) 
   20   continue
      else
        do 30 i=1,nvol              
          const2(i) = 1.d0
   30   continue
      endif
c
      do 40 i=1,nvol
         q(i,1) = const2(i)*qaux(i,1) - const*rhs(i,1)
         q(i,2) = const2(i)*qaux(i,2) - const*rhs(i,2)
         q(i,3) = const2(i)*qaux(i,3) - const*rhs(i,3)
         q(i,4) = const2(i)*qaux(i,4) - const*rhs(i,4)
   40 continue
c
c     Now, do Stages # 2 to 5.
c
      do 200 k=2,5
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
c        In the 2nd stage, we must compute both convective and
c        dissipative operators. In the other stages (i.e., 3-5),
c        only the convective operator must be computed.
c
         if(k.eq.2) then
            call convec
            if (idissip.eq.0) then
               call dissip
            elseif (idissip.eq.1) then
               call dissmavri
            endif
         else
            call convec
         end if
c
c        Perform the Runge-Kutta stage.
c
         const = h * alphp(k)
c
         if (iflow2.eq.1 .and. niter.gt.1) then
           do 50 i=1,nvol
             const1(i) = const/vol(i,2)                  
   50      continue
         else
           do 60 i=1,nvol              
             const1(i) = const/vol(i,1)
   60      continue
         endif
c
         do 150 i=1,nvol
            c1 = const2(i)*qaux(i,1)
            c2 = const2(i)*qaux(i,2)
            c3 = const2(i)*qaux(i,3)
            c4 = const2(i)*qaux(i,4)
c 
            q(i,1) = c1 - const1(i)*(cpen(i,1)-dpen(i,1))
            q(i,2) = c2 - const1(i)*(cpen(i,2)-dpen(i,2))
            q(i,3) = c3 - const1(i)*(cpen(i,3)-dpen(i,3))
            q(i,4) = c4 - const1(i)*(cpen(i,4)-dpen(i,4))
  150    continue
  200 continue
c
      return
c
c     This is the case with implicit residual smoothing (or implicit
c     residual averaging) but without variable time step.
c
  300 continue
c
c     The quantity that comes into this routine is the current
c     unsmoothed residue. Hence, we must smooth it first.
c
      call resmth
c
c     Stage # 1.
c
      const = h * alphp(1)
      do 350 i=1,nvol
         q(i,1) = qaux(i,1) - const * rhs(i,1)
         q(i,2) = qaux(i,2) - const * rhs(i,2)
         q(i,3) = qaux(i,3) - const * rhs(i,3)
         q(i,4) = qaux(i,4) - const * rhs(i,4)
  350 continue
c
c     Now, do Stages # 2 to 5.
c
      do 400 k=2,5
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
c        In the 2nd stage, we must compute both convective and
c        dissipative operators. In the other stages (i.e., 3-5),
c        only the convective operator must be computed.
c
         if(k.eq.2) then
            call convec
            if (idissip.eq.0) then
               call dissip
            elseif (idissip.eq.1) then
               call dissmavri
            endif
         else
            call convec
         end if
c
c        Compute the new residue (unsmoothed).
c
         do 370 i=1,nvol
            rhs(i,1) = (cpen(i,1)-dpen(i,1)) / vol(i,1)
            rhs(i,2) = (cpen(i,2)-dpen(i,2)) / vol(i,1)
            rhs(i,3) = (cpen(i,3)-dpen(i,3)) / vol(i,1)
            rhs(i,4) = (cpen(i,4)-dpen(i,4)) / vol(i,1)
  370    continue
c
c        Smooth it when appropriate.
c
         if(k.eq.3 .or. k.eq.5) then
            call resmth
         else
            if(irstag.eq.1) call resmth
         end if
c
c        Perform the Runge-Kutta stage.
c
         const = h * alphp(k)
         do 390 i=1,nvol
            q(i,1) = qaux(i,1) - const * rhs(i,1)
            q(i,2) = qaux(i,2) - const * rhs(i,2)
            q(i,3) = qaux(i,3) - const * rhs(i,3)
            q(i,4) = qaux(i,4) - const * rhs(i,4)
  390    continue
  400 continue
c
      return
c
c     This is the variable time step case.
c
c     If this option is selected, remember that "h" is the CFL number.
c
  500 continue
c
c     Here, we must also decide whether there is going to be residual
c     smoothing or not.
c
      if(irsmth.eq.1) go to 800
c
c     OK, then it is going to be variable time stepping but without
c     residual smoothing.
c
c     The quantity that comes into this routine is already the current
c     residue. Therefore, we can go straight to stage # 1.
c
c     Stage # 1.
c
      do 550 i=1,nvol
         const = alphp(1) * deltat(i)
         q(i,1) = qaux(i,1) - const*rhs(i,1)
         q(i,2) = qaux(i,2) - const*rhs(i,2)
         q(i,3) = qaux(i,3) - const*rhs(i,3)
         q(i,4) = qaux(i,4) - const*rhs(i,4)
  550 continue
c
c     Now, do Stages # 2 to 5.
c
      do 700 k=2,5
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
c        In the 2nd stage, we must compute both convective and
c        dissipative operators. In the other stages (i.e., 3-5),
c        only the convective operator must be computed.
c
         if(k.eq.2) then
            call convec
            if (idissip.eq.0) then
               call dissip
            elseif (idissip.eq.1) then
               call dissmavri
            endif
         else
            call convec
         end if
c
c        Perform the Runge-Kutta stage.
c
         do 650 i=1,nvol
            const = alphp(k) * deltat(i) / vol(i,1)
            q(i,1) = qaux(i,1) - const*(cpen(i,1)-dpen(i,1))
            q(i,2) = qaux(i,2) - const*(cpen(i,2)-dpen(i,2))
            q(i,3) = qaux(i,3) - const*(cpen(i,3)-dpen(i,3))
            q(i,4) = qaux(i,4) - const*(cpen(i,4)-dpen(i,4))
  650    continue
  700 continue
c
      return
c
c     Finally, this is the case with variable time stepping and with
c     implicit residual smoothing (or implicit residual averaging).
c
  800 continue
c
c     The quantity that comes into this routine is the current
c     unsmoothed residue. Hence, the first step is to smooth it,
c     before performing the first stage of the R-K time stepping.
c
      call resmth
c
c     Stage # 1.
c
      do 850 i=1,nvol
         const = alphp(1) * deltat(i)
         q(i,1) = qaux(i,1) - const*rhs(i,1)
         q(i,2) = qaux(i,2) - const*rhs(i,2)
         q(i,3) = qaux(i,3) - const*rhs(i,3)
         q(i,4) = qaux(i,4) - const*rhs(i,4)
  850 continue
c
c     Now, perform Stages # 2 to 5.
c
      do 1000 k=2,5
c
c        It is necessary to update the b.c.'s during the Runge-Kutta
c        stages. Hence, call subroutine Boundary.
c
         call boundary
c
c        In the 2nd stage, we must compute both convective and
c        dissipative operators. In the other stages (i.e., 3-5),
c        only the convective operator must be computed.
c
         if(k.eq.2) then
            call convec
            if (idissip.eq.0) then
               call dissip
            elseif (idissip.eq.1) then
               call dissmavri
            endif
         else
            call convec
         end if
c
c        Compute the new (unsmoothed) residue using the current fluxes.
c
         do 900 i=1,nvol
            rhs(i,1) = (cpen(i,1) - dpen(i,1)) / vol(i,1)
            rhs(i,2) = (cpen(i,2) - dpen(i,2)) / vol(i,1)
            rhs(i,3) = (cpen(i,3) - dpen(i,3)) / vol(i,1)
            rhs(i,4) = (cpen(i,4) - dpen(i,4)) / vol(i,1)
  900    continue
c
c        Perform the appropriate residual smoothing for this stage.
c
         if(k.eq.3 .or. k.eq.5) then
            call resmth
         else
            if(irstag.eq.1) call resmth
         end if
c
c        Perform the Runge-Kutta stage.
c
         do 950 i=1,nvol
            const = alphp(k) * deltat(i)
            q(i,1) = qaux(i,1) - const * rhs(i,1)
            q(i,2) = qaux(i,2) - const * rhs(i,2)
            q(i,3) = qaux(i,3) - const * rhs(i,3)
            q(i,4) = qaux(i,4) - const * rhs(i,4)
  950    continue
 1000 continue
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine START
c
c     This subroutine initializes the flow variables for every control
c     volume. If this is a restart case, the flow variables are simply
c     read from a file which contains the results of the previous run
c     (remember that, in this case, file Fort.7 must exist).
c
      subroutine start
c
      include 'Common_rm_v2.0.f'
c
      if(iflow2.eq.0) then
c
c     This is the steady case. First of all, we must zero out the
c     components of the grid and cell velocities.
c            
         alphat = 0.d0
         ht = 0.d0
c
         do 120 i=1,nsides
            vside(i,1) = 0.d0
            vside(i,2) = 0.d0
  120    continue
c
         do 130 i=1,nvolg
            uvolg(i) = 0.d0
            vvolg(i) = 0.d0
  130    continue      
c
         finaltime = 0.d0
         open(20,file='Fort.20')
c
         if(nstart.eq.0) then
c
c     If this is a true start, generate some initial solution.
c     The initial solution generated depends on the type of case
c     treated in the present run. This is defined by the flag
c     variable "iflow". The possible options implemented so far
c     are:
c               iflow = 1  >  external flow
c                     = 2  >  internal flow, supersonic entrance
c                     = 3  >  internal flow, subsonic entrance
c
c     For cases iflow=1 or =2, freestream conditions are assumed
c     as initial condition (where, for case iflow=2, this means
c     that the entrance condition is assigned to the whole field).
c     For case iflow=3, stagnation conditions are used.
c        
            write(14,100)
c
c     Freestream initial conditions.
c
            if(iflow.le.2) then
               rhoinf = 1.0d0
               uinf = fsmach * dcos(alpha)
               vinf = fsmach * dsin(alpha)
               pinf = 1.0d0 / gamma
               einf = (pinf/(gamma - 1.0d0)) + 0.5d0*rhoinf*(fsmach**2)
c
               aux1 = rhoinf * uinf
               aux2 = rhoinf * vinf
               do 210 i=1,nvol
                  q(i,1) = rhoinf
                  q(i,2) = aux1
                  q(i,3) = aux2
                  q(i,4) = einf
  210          continue
c
c        Stagnation initial conditions.
c
            else
               rhoinf = 1.0d0
               uinf = 0.0d0
               vinf = 0.0d0
               einf = rhoinf * ((gamma + 1.0d0)/
     &                (2.0d0*gamma*(gamma - 1.0d0)) 
     &              + 0.5d0*(uinf*uinf + vinf*vinf))
c
               do 220 i=1,nvol
                  q(i,1) = rhoinf
                  q(i,2) = 0.0d0
                  q(i,3) = 0.0d0
                  q(i,4) = einf
  220          continue
            endif
         else
c
c     If this is a restart case, read the previous solution from
c     unit = 7. In this case, you also have to open the appropriate
c     data file.
c
            open(7,file='Fort.7',form='unformatted')
            read(7) nstart
            read(7) ((q(i,n),i=1,nvol),n=1,4)
            close(7)
            cltemp = 0.d0
            clsdtemp = 0.d0
            cdtemp = 0.d0
            cmtemp = 0.d0
            cmsdtemp = 0.d0
            cmextemp = 0.d0
            read(14,*)
            naux = nstart/nhist
            do i = 1,naux
               read(14,900) temp1,temp2,temp3,dcl,dcd,dcm,dclsd,dcmsd,
     &                      dcmex
               cltemp = cltemp + dcl
               clsdtemp = clsdtemp + dclsd
               cdtemp = cdtemp + dcd
               cmtemp = cmtemp + dcm
               cmsdtemp = cmsdtemp + dcmsd
               cmextemp = cmextemp + dcmex
            enddo
            do i = 1,nstart 
               read(9,*)
               read(11,*)
            enddo
            open(22,file='Fort.22')
            read(22,*)
            read(22,950) finaltime
            close(22)
         endif
c
c     The boundary conditions on the initial data must be imposed at
c     this point. Observe that we are assuming that b.c.'s will be
c     imposed with the help of "ghost" control volumes (actually,
c     MacCormack calls then "slave" control volumes). Moreover, the
c     information in these volumes is not saved from one run to the
c     other. Therefore, it must be reconstructed using the interior
c     information.
c
            call boundary
c
      else
c
c     This is the unsteady case. Now we first set the number of
c     iterations where mesh motion is required.
c
         if (imove.eq.0) then
            nref = nmiter + 1
         elseif (imove.eq.1) then
            nref = tmax/h
            nref = nref + 1
         elseif (imove.eq.2) then
            nref = 2
         else
            if (imotion.lt.2) then
               nref = 3
            else
               nref = 2
            endif
         endif
c
c     If intermediary results are to be saved, prepare Fort.23
c    
         if(iresp.eq.1) then
            open(23,file='Fort.23',form='unformatted')
            write(23) nnodes,gamma,fsmach,iflow
         endif
         if(ipress.eq.1) open(19,file='Fort.19')
c       
         if(nstart.eq.0) then
c
c     This is an unsteady start case, read the the steady-state
c     solution from  unit = 7. In this case, you also have to open
c     the appropriate data file. It's important to say that we still
c     perform another 'steady' iteration before going to move the 
c     airfoil. We must also zero out the components of the grid and
c     cell velocities.
c            
            alphat = 0.d0
            ht = 0.d0
c
            do 320 i =1,nsides
               vside(i,1) = 0.d0
               vside(i,2) = 0.d0
  320       continue
            do 330 i=1,nvolg
               uvolg(i) = 0.d0
               vvolg(i) = 0.d0
  330       continue      
c
            cltemp = 0.d0
            clsdtemp = 0.d0
            cdtemp = 0.d0
            cmtemp = 0.d0
            cmsdtemp = 0.d0
            cmextemp = 0.d0
            write(14,100)
            finaltime = 0.d0
            open(7,file='Fort.7',form='unformatted')
            read(7) ndum
            read(7) ((q(i,n),i=1,nvol),n=1,4)
            close(7)
            open(21,file='Fort.21')
            call boundary
c
         else
c
c     This is the unsteady flow restart implementation 
c
            open(15,file='Fort.15',form='unformatted')
c
            read(15) nstart
            read(15) ((q(i,n),i=1,nvol),n=1,4)
            read(15) ((rghost(i,n),i=1,nvolg),n=1,5)
            read(15) ((vol(i,n),i=1,nvol),n=1,2)
            read(15) ((xynode(i,n),i=1,nnodes),n=1,2)
            read(15) ((xhist(i,n),i=1,nnodes),n=1,2)
            read(15) ((yhist(i,n),i=1,nnodes),n=1,2)
            read(15) ((vside(i,n),i=1,nsides),n=1,2)
            read(15) (uvolg(i),i=1,nvolg)
            read(15) (vvolg(i),i=1,nvolg)
            read(15) alphat, ht
            read(15) (x0(i),i=1,nnodes)
            read(15) (y0(i),i=1,nnodes)
            close(15)
c
            cltemp = 0.d0
            clsdtemp = 0.d0
            cdtemp = 0.d0
            cmtemp = 0.d0
            cmsdtemp = 0.d0
            cmextemp = 0.d0
            read(14,*)
            naux = nstart/nhist
            do i = 1,naux
               read(14,900) temp1,temp2,temp3,dcl,dcd,dcm,dclsd,dcmsd,
     &                      dcmex
               cltemp = cltemp + dcl
               clsdtemp = clsdtemp + dclsd
               cdtemp = cdtemp + dcd
               cmtemp = cmtemp + dcm
               cmsdtemp = cmsdtemp + dcmsd
               cmextemp = cmextemp + dcmex
            enddo
            do i = 1,nstart 
               read(9,*)
               read(11,*)
            enddo
            open(22,file='Fort.22')
            read(22,*)
            read(22,950) finaltime
            close(22)
            if(nstart.lt.(nref-1)) then
               open(21,file='Fort.21')
               do i = 1,nstart
                  read(21,*)
               enddo
            endif              
        endif
      endif
c
      return
c
c     Formats.
c
100   format('time ht alphai cl cd cm clsd cmsd cmex')
900   format(9(e12.6,1x))
950   format(62x,e9.3)
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine STOPIT
c
c     This subroutine stops the code according to the value of icode.
c     If:  icode = 0 > means that the maximum number of iterations has
c                      been exceeded; the code is stopped and the
c                      current solution is saved.
c                = 1 > means that the code has converged; stop and
c                      save final result.
c                = 2 > the code may be diverging; it is stopped, but
c                      the current solution is not saved (idea is that
c                      one could use the previously saved solution,
c                      before blowing up, to restart).
c                = 3 > actually, it does not stop the code in this
c                      case, but it simply saves the current solution
c                      in the appropriate file and returns to main.
c
      subroutine stopit(icode)
c
      include 'Common_rm_v2.0.f'
c
c     Define some temporary storage.
c
      real*4 addtime,timecomp(2)
c      common p(2,nmaxgh),x(2,nmaxgh),y(2,nmaxgh),
c     &          dens(2,nmaxgh),entr(2,nmaxgh)
c      common iup,ilow
c
c     Unless the code has diverged or we are only onterested in pressure
c     distribution, the solution has to be saved.
c
      if (icode.ne.2) then
c
c     For the steady case, save the current solution in a format ready
c     for input for the next run.
c
         if(iflow2.eq.0) then
            open(8,file='Fort.8',form='unformatted')
            write(8) niter
            write(8) ((q(i,n),i=1,nvol),n=1,4)
            write(8) gamma,fsmach,iflow
            close(8)
            addtime = dtime(timecomp)
            finaltime = addtime + finaltime
            open(22,file='Fort.22')
            read(22,*)
            write(22,930) finaltime
            close(22)
            close(9)
            close(11)
            close(14)
            if (icode.eq.3) then
               open(9,file='Fort.9')
               open(11,file='Fort.11')
               open(14,file='Fort.14')
               do i = 1,niter
                  read(9,*)
                  read(11,*)
               enddo
               naux = niter/nhist
               do i = 1,(naux+1)
                  read(14,*)
               enddo
            endif
         else
c
c      These are the unsteady flow restart saves  
c
            open(16,file='Fort.16',form='unformatted')
            write(16) niter
            write(16) ((q(i,n),i=1,nvol),n=1,4)
            write(16) ((rghost(i,n),i=1,nvolg),n=1,5)
            write(16) ((vol(i,n),i=1,nvol),n=1,2)
            write(16) ((xynode(i,n),i=1,nnodes),n=1,2)
            write(16) ((xhist(i,n),i=1,nnodes),n=1,2)
            write(16) ((yhist(i,n),i=1,nnodes),n=1,2)
            write(16) ((vside(i,n),i=1,nsides),n=1,2)
            write(16) (uvolg(i),i=1,nvolg)
            write(16) (vvolg(i),i=1,nvolg)
            write(16) alphat, ht
            write(16) (x0(i),i=1,nnodes)
            write(16) (y0(i),i=1,nnodes)
            close(16)         
            addtime = dtime(timecomp)
            finaltime = addtime + finaltime
            open(22,file='Fort.22')
            read(22,*)
            write(22,930) finaltime
            close(22)
            close(9)
            close(11)
            close(14)
            open(9,file='Fort.9')
            open(11,file='Fort.11')
            open(14,file='Fort.14')
            do i = 1,niter
               read(9,*)
               read(11,*)
            enddo
            naux = niter/nhist
            do i = 1,(naux+1)
               read(14,*)
            enddo
            if(niter.lt.(nref-1)) then
               close(21)
               open(21,file='Fort.21')      
               do i = 1,niter
                  read(21,*)
               enddo
            endif
         endif
      endif
c
c     Check if this is a stop case or not.
c
      if(icode .eq. 3) then
         return
      else
         continue
      endif
c
c     Now, before stopping, calculate the wall pressure distribution
c     and the aerodynamic coefficients and print them with the
c     appropriate convergence information.
c
      if(iflow2.eq.0) then
         if(mod(nmiter,nhist).ne.0) call history
         if(mod(nmiter,npressave).ne.0
     &      .or. ipress.eq.0) call pressure
      endif
      if(icode.eq.0) then
        write(*,910) nmiter,rhsmax
      elseif(icode.eq.1) then
        naux = niter - 1
        write(*,921) naux,rhsmax
      elseif(icode.eq.2) then
        write(*,922) niter,rhsmax  
      endif
c
c     Write the elapsed calculation time, discouting
c     the preprocessing phase. 
c
      close(9)
      close(11)
      close(14)
      close(22)
      stop
c 
c     Formats.
c
  910 format(/////,1x,'STOPIT, code has exceeded maximum number of ',
     & 'iterations specified',/,1x,'(Nmax = ',i6,') without achieving',
     & ' convergence (RHSmax = ',e15.7,')')
  921 format(/////,1x,'STOPIT, code converged to required accuracy ',
     & 'in ',i6,' iterations.',/,1x,'RHSmax = ',e15.7)
  922 format(/////,1x,'STOPIT, code seems to be diverging at ',
     & 'iteration number ',i6,/,1x,'RHSmax = ',e15.7) 
  930 format(1x,'Simulation time, discouting the preprocessing phase',
     & ' (CPU s):',1x,e9.3)
c
      end
c ---------------------------------------------------------------------
c
c     Subroutine INTERMED
c
c     This subroutine saves the intermediary results for plotting.
c
      subroutine intermed
c
      include 'Common_rm_v2.0.f'
c      
      time = (niter-1)*h                           
      write(23) time,ht,alphat
      write(23) ((q(i,n),i=1,nvol),n=1,4)
      write(23) ((xynode(i,n),i=1,nnodes),n=1,2)
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine UPDATE
c
c     This subroutine performs the update of all interior volumes and
c     the update of all boundary "ghost" volumes. It also performs the
c     displacement of the mesh nodes for the unsteady flow case.
c
      subroutine update
c
      include 'Common_rm_v2.0.f'
c
c     Actually, there is no need to perform any update of interior
c     volume information. The "q" values that come out from the
c     RK5S subroutine are already the new values at time (n+1).
c     Therefore, in this sense, this subroutine could be eliminated
c     and the main program could call subroutine Boundary directly.
c     However, I decided against eliminating this subroutine at
c     this point (Oct./23/91) because I may decide to implement an
c     implicit solver in the future. In that case, the program
c     would already be prepared for such a modification. If I really
c     stay with the hybrid explicit time-stepping schemes I am
c     using right now, this subroutine should be eliminated in the
c     future.
c
c     The "ghost" volumes are taken care of by subroutine Boundary.
c
c      call boundary
c                                                               
c     For the unsteady flow case, we have to move the mesh before the
c     next iteration. Except for the step case, where the mesh is
c     moved only in the first iteration.
c
      if(iflow2.eq.1) then 
         if(niter.lt.nref) then
c                
c        First we store the volumes (or areas) of each cell.
c        
            if(niter.gt.1) then
               do 40 i=1,nvol
                  vol(i,1) = vol(i,2)
   40          continue
            endif
c
c        And fill the appropriate portion of the "rghost" array with
c        the "volume" (or "area") of that ghost volume.
c
            do 50 ig=1,nvolg
               intngb = ighost(ig,2)
               rghost(ig,5) = vol(intngb,1)
   50       continue
c
c         Then, we move the mesh and calculate the mesh velocities.
c
            if (imotion.eq.0 .or. imotion.eq.1) then
               if (imove.ne.1 .and. imove.ne.3) then
                  call dynmesh
               elseif (niter.lt.(nref-1)) then
                  call dynmesh
               else
                  alphat = 0.d0
                  ht = 0.d0
                  do 55 i=1,nnodes
                     xynode(i,1) = x0(i)
                     xynode(i,2) = y0(i)
   55             continue
               endif
               if (imotion.eq.0) call velocity
            else
               if (niter.gt.1) then
                  do 57 i = 1,nnodes
                     xynode(i,1) = xytemp(i,1)
                     xynode(i,2) = xytemp(i,2)
   57             continue
               endif
               call dynmesh
               call velocity
             endif          	       
c
c     Compares the new mesh with the initial (0)
c     and previous ones
c
            dmesh0m = 0.d0
            dmesh0g = 0.d0
            dmeshm = 0.d0
            dmeshg = 0.d0
            do 160 i=1,nnodes
               dxmesh0 = (xynode(i,1) - x0(i))**2
               dymesh0 = (xynode(i,2) - y0(i))**2
               dmesh0i = dsqrt(dxmesh0 + dymesh0)
               dmesh0g = (dmesh0g + dmesh0i +
     &                    dabs(dmesh0g - dmesh0i))/2.d0
               dmesh0m = dmesh0m + dmesh0i
               dxmesh = (xynode(i,1) - xhist(i,1))**2
               dymesh = (xynode(i,2) - yhist(i,1))**2
               dmeshi = dsqrt(dxmesh + dymesh)
               dmeshg = (dmeshg + dmeshi +
     &                   dabs(dmeshg - dmeshi))/2.d0
               dmeshm = dmeshm + dmeshi
  160       continue
            dmesh0m = dmesh0m/nnodes
            dmeshm = dmeshm/nnodes
            write(21,900) niter,dmesh0m,dmesh0g,dmeshm,dmeshg    
            
            if (imotion.eq.2) then
               do 59 i=1,nnodes
                  xytemp(i,1) = xynode(i,1)
                  xytemp(i,2) = xynode(i,2)
                  xynode(i,1) = x0(i)
                  xynode(i,2) = y0(i)
   59          continue
            endif
c
c         And now, we store the new volumes.
c
            if (igcl.lt.1) then
               do 60 i=1,nvol
                  ina = itable(i,1)
                  inb = itable(i,2)
                  inc = itable(i,3)
c
                  ax = xynode(ina,1)
                  ay = xynode(ina,2)
                  bx = xynode(inb,1)
                  by = xynode(inb,2)
                  cx = xynode(inc,1)
                  cy = xynode(inc,2)
c
                  call volumes(ax,ay,bx,by,cx,cy,area)
                  vol(i,2) = area
   60          continue
            else
               call gcl
            endif
c
            do 70 ig=1,nvolg
               intngb = ighost(ig,2)
               rghost(ig,6) = vol(intngb,2)
   70       continue
c
         elseif(niter.eq.nref) then
            close(21)
            if (imove.ne.2 .or. imotion.ne.2) then
               do 80 i=1,nsides
                  vside(i,1) = 0.d0
                  vside(i,2) = 0.d0
   80          continue
               do 90 i=1,nvol
                  vol(i,1) = vol(i,2)
   90          continue              
               do 100 ig=1,nvolg
                  intngb = ighost(ig,2)
                  rghost(ig,5) = vol(intngb,1)
                  rghost(ig,6) = vol(intngb,2)
                  uvolg(ig) = 0.d0
                  vvolg(ig) = 0.d0
  100          continue
            else
c
c      Case imotion = 2 and imove =2, the velocity still acts on the
c      problem, and, therefore, on the volumes if gcl = 1
c
               do 110 i=1,nvol
                  vol(i,1) = vol(i,2)
  110          continue
c              
               if (igcl.lt.1) then
                  do 120 i=1,nvol
                     vol(i,2) = vol(i,1)
  120             continue
               else
                  call gcl
               endif
               do 130 ig=1,nvolg
                  intngb = ighost(ig,2)
                  rghost(ig,5) = vol(intngb,1)
                  rghost(ig,6) = vol(intngb,2)
  130          continue
            endif
         endif
      endif
c
c     Format
c
  900 format(i6,4(1x,e12.6))
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine GCL (Fred - 05-oct-2003)
c
c     This subroutine computes the new areas for all the "real" volumes
c     using GCL ideas as presented by Batina (AIAA-89-0115).
c
c     Here I am again taking full benefit of the edge based data.
c
      subroutine gcl
c
      include   'Common_rm_v2.0.f'
c
      integer   i, n1, n2, v1, v2
c
      real*8    delx, dely, dxdt, dydt, expr
c
c     Cells loop to initialize new volumes
c
      do i=1,nvol
         vol(i,2) = vol(i,1)
      enddo
c
c     Edges loop to compute sums
c
      do i=1,nsides
         n1 = iedges(i,1)
         n2 = iedges(i,2)
         v1 = iedges(i,3)
         v2 = iedges(i,4)
c
         if (v2.lt.0) then
            delx = xynode(n2,1) - xynode(n1,1)
            dely = xynode(n2,2) - xynode(n1,2)
            expr = h * (vside(i,1)*dely - vside(i,2)*delx)
            vol(v1,2) = vol(v1,2) + expr
         else
            delx = xynode(n2,1) - xynode(n1,1)
            dely = xynode(n2,2) - xynode(n1,2)
            expr = h * (vside(i,1)*dely - vside(i,2)*delx)
            vol(v1,2) = vol(v1,2) + expr
            vol(v2,2) = vol(v2,2) - expr
         endif
      enddo
c
      return
      end
c ---------------------------------------------------------------------
c
c     Subroutine VARTIME
c
c     This subroutine computes the space-varying time step. This
c     calculation is performed using the characteristic length for
c     each computational cell, computed in Subroutine Meshio, and the
c     current value of the maximum possible signal propagation speed.
c     This is assumed to be (u_total + aa), where u_total is the
c     magnitude of velocity vector computed with the cell properties
c     and aa is the local speed of sound. The CFL number is assumed
c     constant, and it must be emphasized that, in this case, the
c     variable "h" is the CFL number.
c
      subroutine vartime
c
      include 'Common_rm_v2.0.f'
c
c     This is mostly a loop over all the computational cells.
c
      const = gamma - 1.0d0
      do 100 i=1,nvol
c
c        Compute pressure.
c
         p = const * (q(i,4) - 0.5d0 * (q(i,2)*q(i,2) +
     &                                  q(i,3)*q(i,3)) / q(i,1))
c
c        Compute speed of sound.
c
         aa = dsqrt(gamma * p / q(i,1))
c
c        Compute magnitude of local velocity vector.
c
         uu = q(i,2) / q(i,1)
         vv = q(i,3) / q(i,1)
         utotal = dsqrt(uu*uu + vv*vv)
c
c        Finally, compute local time step.
c
         deltat(i) = h * clength(i) / (utotal + aa)
  100 continue
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine VLFLUX
c
c     This subroutine computes the E_plus, E_minus, F_plus and F_minus
c     flux vectors, according to the Van Leer flux vector splitting
c     scheme, using the properties of each control cell. I emphasize,
c     the calculation is performed in a control volume basis, that is,
c     the upwind fluxes are calculated for each control volume using
c     (obviously) the properties of that control volume. This is not
c     an interface flux calculation.
c
c     The computation is performed for both "interior" (or "real") and
c     "ghost" volumes.
c
      subroutine vlflux
c
      include 'Common_rm_v2.0.f'
c
c     Some temporary storage.
c
      dimension u(nmaxvl),v(nmaxvl),p(nmaxvl),a(nmaxvl)
c
c     First, perform the calculation for the "real" volumes.
c
c     Perform preliminary computations.
c
      const = gamma - 1.0d0
      const2 = (gamma * gamma) - 1.0d0
      do 20 i=1,nvol
         u(i) = q(i,2) / q(i,1)
         v(i) = q(i,3) / q(i,1)
         p(i) = const * (q(i,4) - 0.5d0*(q(i,2)*q(i,2)
     &                                 + q(i,3)*q(i,3))/q(i,1))
   20 continue
c
      do 30 i=1,nvol
         a(i) = dsqrt((gamma * p(i))/q(i,1))
   30 continue
c
c     Compute the E and F flux vectors.
c
      do 100 i=1,nvol
c
c        Do the E_plus and E_minus vectors first.
c
         amx = u(i) / a(i)
         if(amx.ge.1.0d0) then
c
c           For Mx > 1: E_plus = E and E_minus = 0.
c
            eplus(i,1) = q(i,2)
            eplus(i,2) = q(i,2)*u(i) + p(i)
            eplus(i,3) = q(i,2)*v(i)
            eplus(i,4) = (q(i,4) + p(i)) * u(i)
            eminus(i,1) = 0.0d0
            eminus(i,2) = 0.0d0
            eminus(i,3) = 0.0d0
            eminus(i,4) = 0.0d0
         else
c
c           For Mx < -1: E_plus = 0 and E_minus = E.
c
            if(amx.le.-1.0d0) then
               eplus(i,1) = 0.0d0
               eplus(i,2) = 0.0d0
               eplus(i,3) = 0.0d0
               eplus(i,4) = 0.0d0
               eminus(i,1) = q(i,2)
               eminus(i,2) = q(i,2)*u(i) + p(i)
               eminus(i,3) = q(i,2)*v(i)
               eminus(i,4) = (q(i,4) + p(i)) * u(i)
c
            else
c
c           For -1 < Mx < 1: use the Van Leer formulas.
c
               fpmass = q(i,1) * a(i) * (0.5d0 * (amx + 1.0d0))
     &                                * (0.5d0 * (amx + 1.0d0))
               fmmass = - q(i,1) * a(i) * (0.5d0 * (amx - 1.0d0))
     &                                  * (0.5d0 * (amx - 1.0d0))
               aux1 = const * u(i) + 2.0d0 * a(i)
               aux2 = const * u(i) - 2.0d0 * a(i)
               eplus(i,1) = fpmass
               eplus(i,2) = fpmass * aux1 / gamma
               eplus(i,3) = fpmass * v(i)
               eplus(i,4) = fpmass * (((aux1*aux1)/(2.0d0*const2))
     &                             + (0.5d0 * v(i) * v(i)))
               eminus(i,1) = fmmass
               eminus(i,2) = fmmass * aux2 / gamma
               eminus(i,3) = fmmass * v(i)
               eminus(i,4) = fmmass * (((aux2*aux2)/(2.0d0*const2))
     &                              + (0.5d0 * v(i) * v(i)))
            end if
         end if
c
c        Now, do the F_plus and F_minus vectors.
c
         amy = v(i) / a(i)
         if(amy.ge.1.0d0) then
c
c           For My > 1: F_plus = F and F_minus = 0.
c
            fplus(i,1) = q(i,3)
            fplus(i,2) = q(i,3)*u(i)
            fplus(i,3) = q(i,3)*v(i) + p(i)
            fplus(i,4) = (q(i,4) + p(i)) * v(i)
            fminus(i,1) = 0.0d0
            fminus(i,2) = 0.0d0
            fminus(i,3) = 0.0d0
            fminus(i,4) = 0.0d0
         else
c
c           For My < -1: F_plus = 0 and F_minus = F.
c
            if(amy.le.-1.0d0) then
               fplus(i,1) = 0.0d0
               fplus(i,2) = 0.0d0
               fplus(i,3) = 0.0d0
               fplus(i,4) = 0.0d0
               fminus(i,1) = q(i,3)
               fminus(i,2) = q(i,3)*u(i)
               fminus(i,3) = q(i,3)*v(i) + p(i)
               fminus(i,4) = (q(i,4) + p(i)) * v(i)
c
            else
c
c           For -1 < My < 1: use the Van Leer formulas.
c
               fpmass = q(i,1) * a(i) * (0.5d0 * (amy + 1.0d0))
     &                                * (0.5d0 * (amy + 1.0d0))
               fmmass = - q(i,1) * a(i) * (0.5d0 * (amy - 1.0d0))
     &                                  * (0.5d0 * (amy - 1.0d0))
               aux1 = const * v(i) + 2.0d0 * a(i)
               aux2 = const * v(i) - 2.0d0 * a(i)
               fplus(i,1) = fpmass
               fplus(i,2) = fpmass * u(i)
               fplus(i,3) = fpmass * aux1 / gamma
               fplus(i,4) = fpmass * (((aux1*aux1)/(2.0d0*const2))
     &                             + (0.5d0 * u(i) * u(i)))
               fminus(i,1) = fmmass
               fminus(i,2) = fmmass * u(i)
               fminus(i,3) = fmmass * aux2 / gamma
               fminus(i,4) = fmmass * (((aux2*aux2)/(2.0d0*const2))
     &                              + (0.5d0 * u(i) * u(i)))
            end if
         end if
c
  100 continue
c
c     Now, perform the calculation for the "ghost" volumes.
c
c     Perform preliminary computations.
c
      const = gamma - 1.0d0
      const2 = (gamma * gamma) - 1.0d0
      do 220 i=1,nvolg
         u(i) = rghost(i,2) / rghost(i,1)
         v(i) = rghost(i,3) / rghost(i,1)
         p(i) = const * (rghost(i,4) - 0.5d0*(rghost(i,2)*rghost(i,2)
     &                 + rghost(i,3)*rghost(i,3))/rghost(i,1))
  220 continue
c
      do 230 i=1,nvolg
         a(i) = dsqrt((gamma * p(i))/rghost(i,1))
  230 continue
c
c     Compute the E and F flux vectors.
c
      do 300 i=1,nvolg
c
c        Do the E_plus and E_minus vectors first.
c
         amx = u(i) / a(i)
         if(amx.ge.1.0d0) then
c
c           For Mx > 1: E_plus = E and E_minus = 0.
c
            egplus(i,1) = rghost(i,2)
            egplus(i,2) = rghost(i,2)*u(i) + p(i)
            egplus(i,3) = rghost(i,2)*v(i)
            egplus(i,4) = (rghost(i,4) + p(i)) * u(i)
            egminus(i,1) = 0.0d0
            egminus(i,2) = 0.0d0
            egminus(i,3) = 0.0d0
            egminus(i,4) = 0.0d0
         else
c
c           For Mx < -1: E_plus = 0 and E_minus = E.
c
            if(amx.le.-1.0d0) then
               egplus(i,1) = 0.0d0
               egplus(i,2) = 0.0d0
               egplus(i,3) = 0.0d0
               egplus(i,4) = 0.0d0
               egminus(i,1) = rghost(i,2)
               egminus(i,2) = rghost(i,2)*u(i) + p(i)
               egminus(i,3) = rghost(i,2)*v(i)
               egminus(i,4) = (rghost(i,4) + p(i)) * u(i)
c
            else
c
c           For -1 < Mx < 1: use the Van Leer formulas.
c
               fpmass = rghost(i,1) * a(i) * (0.5d0 * (amx + 1.0d0))
     &                                     * (0.5d0 * (amx + 1.0d0))
               fmmass = - rghost(i,1) * a(i) * (0.5d0 * (amx - 1.0d0))
     &                                       * (0.5d0 * (amx - 1.0d0))
               aux1 = const * u(i) + 2.0d0 * a(i)
               aux2 = const * u(i) - 2.0d0 * a(i)
               egplus(i,1) = fpmass
               egplus(i,2) = fpmass * aux1 / gamma
               egplus(i,3) = fpmass * v(i)
               egplus(i,4) = fpmass * (((aux1*aux1)/(2.0d0*const2))
     &                              + (0.5d0 * v(i) * v(i)))
               egminus(i,1) = fmmass
               egminus(i,2) = fmmass * aux2 / gamma
               egminus(i,3) = fmmass * v(i)
               egminus(i,4) = fmmass * (((aux2*aux2)/(2.0d0*const2))
     &                               + (0.5d0 * v(i) * v(i)))
            end if
         end if
c
c        Now, do the F_plus and F_minus vectors.
c
         amy = v(i) / a(i)
         if(amy.ge.1.0d0) then
c
c           For My > 1: F_plus = F and F_minus = 0.
c
            fgplus(i,1) = rghost(i,3)
            fgplus(i,2) = rghost(i,3)*u(i)
            fgplus(i,3) = rghost(i,3)*v(i) + p(i)
            fgplus(i,4) = (rghost(i,4) + p(i)) * v(i)
            fgminus(i,1) = 0.0d0
            fgminus(i,2) = 0.0d0
            fgminus(i,3) = 0.0d0
            fgminus(i,4) = 0.0d0
         else
c
c           For My < -1: F_plus = 0 and F_minus = F.
c
            if(amy.le.-1.0d0) then
               fgplus(i,1) = 0.0d0
               fgplus(i,2) = 0.0d0
               fgplus(i,3) = 0.0d0
               fgplus(i,4) = 0.0d0
               fgminus(i,1) = rghost(i,3)
               fgminus(i,2) = rghost(i,3)*u(i)
               fgminus(i,3) = rghost(i,3)*v(i) + p(i)
               fgminus(i,4) = (rghost(i,4) + p(i)) * v(i)
c
            else
c
c           For -1 < My < 1: use the Van Leer formulas.
c
               fpmass = rghost(i,1) * a(i) * (0.5d0 * (amy + 1.0d0))
     &                                     * (0.5d0 * (amy + 1.0d0))
               fmmass = - rghost(i,1) * a(i) * (0.5d0 * (amy - 1.0d0))
     &                                       * (0.5d0 * (amy - 1.0d0))
               aux1 = const * v(i) + 2.0d0 * a(i)
               aux2 = const * v(i) - 2.0d0 * a(i)
               fgplus(i,1) = fpmass
               fgplus(i,2) = fpmass * u(i)
               fgplus(i,3) = fpmass * aux1 / gamma
               fgplus(i,4) = fpmass * (((aux1*aux1)/(2.0d0*const2))
     &                              + (0.5d0 * u(i) * u(i)))
               fgminus(i,1) = fmmass
               fgminus(i,2) = fmmass * u(i)
               fgminus(i,3) = fmmass * aux2 / gamma
               fgminus(i,4) = fmmass * (((aux2*aux2)/(2.0d0*const2))
     &                               + (0.5d0 * u(i) * u(i)))
            end if
         end if
c
  300 continue
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine VOLUME
c
c     This subroutine computes the area of a triangle of vertices A, B
c     and C, or, (ax,ay), (bx,by) and (cx,cy).
c
      subroutine volumes(ax,ay,bx,by,cx,cy,area)
c
      implicit real*8 (a-h,o-z)
c
      aux1 = bx - ax
      aux2 = cy - ay
      aux3 = cx - ax
      aux4 = by - ay
      area = 0.5d0 * (aux1*aux2 - aux3*aux4)
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine HISTORY
c
c     This subroutine calculates the aerodynamic coefficients at each
c     iteration and save them into the file Fort.14. This information 
c     will be used to obtain the Fast Fourrier Transform of the  res-
c     ponse.
c
      subroutine history
c
      include 'Common_rm_v2.0.f'
c
      dimension pall(nmaxgh),p(nmaxgh),xm(nmaxgh),dxw(nmaxgh),
     &          dyw(nmaxgh),sl(nmaxgh),sd(nmaxgh),ym(nmaxgh)
c
      cl = 0.d0
      cm = 0.d0
      cd = 0.d0
      clsd = 0.d0
      cmsd = 0.d0
      cmex = 0.d0
      xc4 = -(0.25d0 + xea) * cos(alphat) + xea
c
c     Identify "ghost" volumes which are actually at a wall and
c     calculate the "area" lengths of the wall segments and the
c     dimensionless pressure.
c
      iwall = 0
      const = gamma - 1.0d0
      do 10 ig=1,nvolg
         ibtype = ighost(ig,1)
         in1 = ighost(ig,3)
         in2 = ighost(ig,4)
         if(ibtype.eq.1) then
            iwall = iwall + 1
            p(iwall) = const * (rghost(ig,4) - 0.5d0*(rghost(ig,2)*
     &      rghost(ig,2) + rghost(ig,3)*rghost(ig,3))/rghost(ig,1))
c
            xm(iwall) = 0.5d0 * (xynode(in1,1) + xynode(in2,1))
            ym(iwall) = 0.5d0 * (xynode(in1,2) + xynode(in2,2))
c
            dxw(iwall) = xynode(in2,1) - xynode(in1,1)
            dyw(iwall) = xynode(in2,2) - xynode(in1,2)
c
            sl(iwall) = dxw(iwall)*cos(alpha) + dyw(iwall)*sin(alpha)
            sd(iwall) = dyw(iwall)*cos(alpha) - dxw(iwall)*sin(alpha)
         end if
   10 continue
c
c     Compute the aerodynamic coefficients.
c
      const = 2.0d0 / (gamma * fsmach * fsmach)
      do 20 i=1,iwall
         p(i) = const * (gamma*p(i) - 1.0d0)
c               
         cl   = cl   - p(i)*sl(i)
         clsd = clsd - p(i)*dxw(i)
         cd   = cd   + p(i)*sd(i)
         cm   = cm   + p(i)*sl(i)*(xm(i)-xc4)
         cmsd = cmsd - p(i)*dxw(i)*(xea-xm(i))
         cmex = cmex - p(i)*dxw(i)*(xea-xm(i))-p(i)*dyw(i)*ym(i)
   20 continue
c
c     Writing into the output file. At this point the mesh
c     motion relative to the present iteration has not been
c     performed and, therefore, the time of the last iteration
c     is accounted.
c     The new way of writting the aerodynamic coefficients is
c     an incremental one. This way, one can get a higher sensibility
c     in unsteady data.
c              
      dcl = cl - cltemp
      dclsd = clsd - clsdtemp
      dcd = cd - cdtemp
      dcm = cm - cmtemp
      dcmsd = cmsd - cmsdtemp
      dcmex = cmex - cmextemp
c
      cltemp = cl
      clsdtemp = clsd
      cdtemp = cd
      cmtemp = cm
      cmsdtemp = cmsd
      cmextemp = cmex
c
      pi = acos(-1.d0)
      alphai = alphat*180.d0/pi
      time = (niter-1)*h                           
c     
      write(14,900) time,ht,alphai,dcl,dcd,dcm,dclsd,dcmsd,
     &              dcmex
c
  900 format(9(e12.6,1x))
      return
c
      end 
c
c ---------------------------------------------------------------------
c
c     Subroutine PREPROC2
c
c     This subroutine identifies the triangles which are to be used
c     in the reconstruction process in order to obtain 2nd order
c     accuracy (in principle, with both the Van Leer and the Liou-
c     Steffen schemes).
c
c     The routine assumes the edge-based data structured has already
c     been set up.
c
      subroutine preproc2
c
      include 'Common_rm_v2.0.f'
c
      dimension itatest(50),itsnode(50)
c
c     This has to be a loop over all edges in order to obtain the
c     information for each of them.
c
      do 500 j=1,nsides
c
c        Get the information from the edge-based database array.
c
         in1 = iedges(j,1)
         in2 = iedges(j,2)
         i = iedges(j,3)
         nb = iedges(j,4)
c
c        First, perform the calculation for the left-hand volume.
c        Remember, this is always a "real" volume by construction.
c
c        Find the centroid.
c
c        (In my notation here, the centroid is called point "n3".)
c        (Moreover, "iotnode" is the other node of triangle "i", i.e.,
c        it is the other node besides "n1" and "n2", which are the
c        two nodes along the edge "j" being currently considered.)
c
         if(in1.eq.itable(i,1)) iotnode = itable(i,3)
         if(in1.eq.itable(i,2)) iotnode = itable(i,1)
         if(in1.eq.itable(i,3)) iotnode = itable(i,2)
c
         xxn3 = (xynode(in1,1)+xynode(in2,1)+xynode(iotnode,1)) / 3.0d0
         yyn3 = (xynode(in1,2)+xynode(in2,2)+xynode(iotnode,2)) / 3.0d0
c
c        Find the point at mid-side.
c        (The mid-side point is "n4".)
c
         xxn4 = 0.5d0 * (xynode(in1,1) + xynode(in2,1))
         yyn4 = 0.5d0 * (xynode(in1,2) + xynode(in2,2))
c
c        Define length of n3-n4 segment and appropriate direction
c        cosines.
c
         aux1 = xxn3 - xxn4
         aux2 = yyn3 - yyn4
         aln3n4 = sqrt((aux1**2) + (aux2**2))
         sx34 = aux1 / aln3n4
         sy34 = aux2 / aln3n4
c
c        Find position of point n5.
c        (More notation, "n5" is the point that is actually going to be
c        used for reconstruction for the left-hand volume.)
c
         aln3n5 = 2.0d0 * aln3n4
         xxn5 = xxn3 + (aln3n5 * sx34)
         yyn5 = yyn3 + (aln3n5 * sy34)
c
c        Now comes the hard part...
c        Identify the triangle in which the point n5 is located.
c
c        First, one has to identify all the triangles that share the
c        node "iotnode".
c
         call findtrn(iotnode,in1,in2,i,itsnode)
c
c        Then, loop over the list of triangles which share "iotnode"
c        (and their neighbors).
c        Remember, the triangles in the "itsnode" list are always
c        "real" triangles. However, their neighbors could be "ghost"
c        triangles.
c
         nctest = 1
         icounter = itsnode(1)
         if(icounter.eq.0) go to 260
c
         do 250 k=2,(icounter+1)
            nkn1l = itsnode(k)
c
c           Perform the area test for that triangle.
c
c           The convention I'll use is:
c               locflag = 1  >  point inside triangle
c               locflag = 0  >  point outside triangle
c
            locflag = 1
            call testinp(nkn1l,locflag,xxn5,yyn5)
c
c           If you found the triangle, write the information in the
c           appropriate array and go to the right-hand volume.
c
            if(locflag.eq.1) then
               iedaux(j,1) = nkn1l
               go to 300
            end if
c
c           If you did not find the triangle, then search its
c           neighbors (remember, these are going to be the neighbors
c           of neighbors).
c
            do 200 kk2=1,3
               nkn2l = neighbor(nkn1l,kk2)
c
c              Eliminate the obvious possibilities.
c
               if(nkn2l.eq.nb .or. nkn2l.eq.i) go to 190
c
c              If the neighbor is a ghost, then do not check anything
c              either.
c
               if(nkn2l.lt.0) go to 190
c
c              Eliminate the ones already tested.
c
               do 170 idummy=1,(nctest-1)
                  if(nkn2l.eq.itatest(idummy)) go to 190
  170          continue
c
c              Finally, also eliminate the ones in the "itsnode" list,
c              because they will be tested anyway (or have already been
c              tested).
c
               do 180 idummy=2,(icounter+1)
                  if(nkn2l.eq.itsnode(idummy)) go to 190
  180          continue
c
c              OK, test the triangle.
c
               locflag = 1
               call testinp(nkn2l,locflag,xxn5,yyn5)
c
c              If you found the triangle, write the information in the
c              appropriate array and go to the right-hand volume.
c
               if(locflag.eq.1) then
                  iedaux(j,1) = nkn2l
                  go to 300
               end if
c
c              If not, add the neighbor (of the neighbor) nkn2l to
c              the list of triangles already tested, and continue
c              the search.
c
               itatest(nctest) = nkn2l
               nctest = nctest + 1
c
  190          continue
  200       continue
c
c           If we get to this point, this is because the triangle nkn1l
c           and its neighbors do not contain point n5. Therefore, add
c           nkn1l to the list of triangles already tested and continue
c           the search.
c
            itatest(nctest) = nkn1l
            nctest = nctest + 1
c
  250    continue
c
c        If we get to this point, then the code could not find the
c        point n5 inside any triangle within the two layers of
c        triangles around "i". Be careful, because this could simply
c        mean that you are too close to the boundary...
c        In any case, I will not try to be fancy here. In this
c        situation, I will simply reduce the scheme to 1st order
c        accuracy.
c
  260    continue
         iedaux(j,1) = i
c
c        Now, let's do the calculation for the right-hand volume.
c
  300    continue
c
c        If nb is a ghost, then there is no computation to be performed
c        because we are going to assume 1st order accuracy at the
c        boundaries.
c
         if(nb.lt.0) go to 480
c
c        Otherwise, do the usual stuff.
c
c        Find the centroid.
c
c        (In my notation here, the centroid of triangle "nb" is called
c        point "n6".)
c        (Moreover, "iotnode" is the other node of triangle "nb", i.e.,
c        it is the other node besides "n1" and "n2", which are the
c        two nodes along the edge "j" being currently considered.)
c
         if(in1.eq.itable(nb,1)) iotnode = itable(nb,2)
         if(in1.eq.itable(nb,2)) iotnode = itable(nb,3)
         if(in1.eq.itable(nb,3)) iotnode = itable(nb,1)
c
         xxn6 = (xynode(in1,1)+xynode(in2,1)+xynode(iotnode,1)) / 3.0d0
         yyn6 = (xynode(in1,2)+xynode(in2,2)+xynode(iotnode,2)) / 3.0d0
c
c        There is no need to find the point at mid-side, because it
c        is the same (obviously)...
c        Remember, mid-side point is "n4".
c
c        Define length of n6-n4 segment and appropriate direction
c        cosines.
c
         aux1 = xxn6 - xxn4
         aux2 = yyn6 - yyn4
         aln6n4 = sqrt((aux1**2) + (aux2**2))
         sx64 = aux1 / aln6n4
         sy64 = aux2 / aln6n4
c
c        Find position of point n7.
c        (More notation, "n7" is the point that is actually going to be
c        used for reconstruction for the right-hand volume.)
c
         aln6n7 = 2.0d0 * aln6n4
         xxn7 = xxn6 + (aln6n7 * sx64)
         yyn7 = yyn6 + (aln6n7 * sy64)
c
c        Now comes the hard part...
c        Identify the triangle in which the point n7 is located.
c
c        As before, one should start by identifying all triangles which
c        share the node "iotnode".
c
         call findtrn(iotnode,in2,in1,nb,itsnode)
c
c        Then, loop over the list of triangles which share "iotnode"
c        (and their neighbors).
c        Remember, the triangles in the "itsnode" list are always
c        "real" triangles. However, their neighbors could be "ghost"
c        triangles.
c
         nctest = 1
         icounter = itsnode(1)
         if(icounter.eq.0) go to 480
c
         do 450 k=2,(icounter+1)
            nkn1l = itsnode(k)
c
c           Perform the area test for that triangle.
c
c           The convention I'll use is:
c               locflag = 1  >  point inside triangle
c               locflag = 0  >  point outside triangle
c
            locflag = 1
            call testinp(nkn1l,locflag,xxn7,yyn7)
c
c           If you found the triangle, write the information in the
c           appropriate array and go to the next edge.
c
            if(locflag.eq.1) then
               iedaux(j,2) = nkn1l
               go to 490
            end if
c
c           If you did not find the triangle, then search its
c           neighbors (remember, these are going to be the neighbors
c           of the neighbor).
c
            do 400 kk2=1,3
               nkn2l = neighbor(nkn1l,kk2)
c
c              Eliminate the obvious possibilities.
c
               if(nkn2l.eq.nb .or. nkn2l.eq.i) go to 390
c
c              If the neighbor is a ghost, then do not check anything
c              either.
c
               if(nkn2l.lt.0) go to 390
c
c              Eliminate the ones already tested.
c
               do 370 idummy=1,(nctest-1)
                  if(nkn2l.eq.itatest(idummy)) go to 390
  370          continue
c
c              Finally, also eliminate the ones in the "itsnode" list,
c              because they will be tested anyway (or have already been
c              tested).
c
               do 380 idummy=2,(icounter+1)
                  if(nkn2l.eq.itsnode(idummy)) go to 390
  380          continue
c
c              OK, test the triangle.
c
               locflag = 1
               call testinp(nkn2l,locflag,xxn7,yyn7)
c
c              If you found the triangle, write the information in the
c              appropriate array and go to the next edge.
c
               if(locflag.eq.1) then
                  iedaux(j,2) = nkn2l
                  go to 490
               end if
c
c              If not, add the neighbor (of the neighbor) nkn2l to
c              the list of triangles already tested, and continue
c              the search.
c
               itatest(nctest) = nkn2l
               nctest = nctest + 1
c
  390          continue
  400       continue
c
c           If we get to this point, this is because the triangle nkn1l
c           and its neighbors do not contain point n7. Therefore, add
c           nkn1l to the list of triangles already tested and continue
c           the search.
c
            itatest(nctest) = nkn1l
            nctest = nctest + 1
c
  450    continue
c
c        If we get to this point, then the code could not find the
c        point n7 inside any triangle within the two layers of
c        triangles around "nb". Be careful, because this could simply
c        mean that you are too close to the boundary...
c        In any case, I will not try to be fancy here. In this
c        situation, I will simply reduce the scheme to 1st order
c        accuracy.
c
  480    continue
         iedaux(j,2) = nb
c
c        OK, we are done with this edge. Go to the next one !!!
c
  490    continue
  500 continue
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine FINDTRN
c
c     This subroutine identifies the triangles which share a given
c     node.
c     One should be careful, however, because this is not a general
c     routine for identifying triangles which share a node. It was
c     coded with the idea of starting at a given triangle, and rotating
c     clockwise around the desired node and marking the triangles in
c     this process. It certainly makes extensive use of (what I call)
c     the neighborhood table (some people call it a "tree").
c     Moreover, if the routine hits a boundary in this process of
c     rotating around the node, it stops and starts at the other edge
c     of the original triangle. If it, then, finds another boundary,
c     it returns for good to the calling routine.
c
c     Input/output data:
c        iinode  > This is the node shared by all the triangles we want
c                  to identify.
c        iednode > This is the node that defines the edge of triangle
c                  "ioritr" where we are going to start the search.
c        i2node  > This is the third node of triangle "ioritr" (needed
c                  here to define the target edge).
c        ioritr  > Identifies the triangle where we are starting the
c                  search.
c        itsnode > This is the desired list. First element of the array
c                  defines how many triangles share node "iinode".
c                  Triangle "ioritr" is NOT included in this list and,
c                  therefore, it is not counted either in "itsnode(1)".
c
      subroutine findtrn(iinode,iednode,i2node,ioritr,itsnode)
c
      include 'Common_rm_v2.0.f'
c
c     Some temporary storage.
c
      dimension itsnode(50)
c
c     Initialize necessary variables.
c
      icount = 0
      icenter = iinode
      iother = iednode
      icurtr = ioritr
c
c     Define the target edge.
c
      itarget = i2node
c
c     Find the next triangle.
c
      if(icenter.eq.itable(icurtr,1)) then
         inextr = neighbor(icurtr,1)
      else
         if(icenter.eq.itable(icurtr,2)) then
            inextr = neighbor(icurtr,2)
         else
            inextr = neighbor(icurtr,3)
         end if
      end if
c
c     Check if next triangle is a ghost.
c
   10 continue
      if(inextr.lt.0) go to 200
c
c     OK, if we get to this point, then it is a real triangle.
c
c     Add triangle to the list of triangles which share "iinode"
c     (or "icenter", if you prefer...).
c
      icount = icount + 1
      itsnode(icount+1) = inextr
c
c     Update variable "icurtr".
c
      icurtr = inextr
c
c     Identify the edge of the new triangle which should be used to
c     continue the search and, since it is the same test, go ahead and
c     already identify the new next triangle (even though you may not
c     need to use this information).
c
      if(icenter.eq.itable(icurtr,1)) then
         iother = itable(icurtr,2)
         inextr = neighbor(icurtr,1)
      else
         if(icenter.eq.itable(icurtr,2)) then
            iother = itable(icurtr,3)
            inextr = neighbor(icurtr,2)
         else
            iother = itable(icurtr,1)
            inextr = neighbor(icurtr,3)
         end if
      end if
c
c     Finally, check whether we are already done with the search.
c
      if(iother.ne.itarget) go to 10
c
c     If we get to this point, it is because your search is over.
c     Write appropriate information in the "itsnode" array and return.
c     (Remember, "icount" is the actual number of triangles you
c     identified and included in "itsnode".)
c
      itsnode(1) = icount
      return
c
c     Now comes the case in which we have hited a boundary before
c     getting to our target edge.
c
  200 continue
c
c     In this case, we have to reverse the direction of the search.
c     Be careful, because things are sort of backwards now...
c     So, start by redefining your control variables.
c
      iother = i2node
      icurtr = ioritr
      itarget = iednode
c
c     Find the next triangle.
c
      if(icenter.eq.itable(icurtr,1)) then
         inextr = neighbor(icurtr,3)
      else
         if(icenter.eq.itable(icurtr,2)) then
            inextr = neighbor(icurtr,1)
         else
            inextr = neighbor(icurtr,2)
         end if
      end if
c
c     Check if next triangle is a ghost.
c
  250 continue
      if(inextr.lt.0) go to 400
c
c     OK, if we get to this point, then it is a real triangle.
c
c     Add triangle to the list of triangles which share "iinode"
c     (or "icenter", if you prefer...).
c
      icount = icount + 1
      itsnode(icount+1) = inextr
c
c     Update variable "icurtr".
c
      icurtr = inextr
c
c     Identify the edge of the new triangle which should be used to
c     continue the search and, since it is the same test, go ahead and
c     already identify the new next triangle (even though you may not
c     need to use this information).
c
      if(icenter.eq.itable(icurtr,1)) then
         iother = itable(icurtr,3)
         inextr = neighbor(icurtr,3)
      else
         if(icenter.eq.itable(icurtr,2)) then
            iother = itable(icurtr,1)
            inextr = neighbor(icurtr,1)
         else
            iother = itable(icurtr,2)
            inextr = neighbor(icurtr,2)
         end if
      end if
c
c     Finally, check whether we are already done with the search.
c
      if(iother.ne.itarget) go to 250
c
c     If we get to this point, it is because your search is over.
c     Write appropriate information in the "itsnode" array and return.
c     (Remember, "icount" is the actual number of triangles you
c     identified and included in "itsnode".)
c     Further, observe that, if we find another boundary, we also
c     return in this case.
c
  400 continue
c
      itsnode(1) = icount
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine LIMITER
c
c     This subroutine computes the limiters that should be used at a
c     given interface (and direction). The implementation here is a
c     straightforward extension of the 1-D ideas to the present
c     2-D, unstructured, triangular case. The previous comment about
c     "direction" has to do with the fact that the routine will have
c     to be called twice, i.e., once when computing the x-direction
c     fluxes and again when computing the y-direction fluxes.
c
c     The limiting is being performed over primitive variables. This
c     is probably not the ideal situation. I would believe that
c     doing the limiting over characteristic variables would be better.
c     However, it would also be more expensive. So, I will do it over
c     the primitive variables...
c
c     The implementation that follows is a preliminary one in the sense
c     that this routine should be generalized in the future in order to
c     allow for the option between a few limiters. In the present
c     version, only the "minmod" limiter is implemented.
c
      subroutine limiter(ilimtyp,qpril1,qpril2,qprir1,qprir2,phimins,
     &                   phiplus)
c
      include 'Common_rm_v2.0.f'
c
      dimension qpril1(4), qpril2(4), qprir1(4), qprir2(4)
      dimension phimins(4), phiplus(4)
      dimension rminus(4), rplus(4)
c
c     Define some preliminary variables. Some of these definitions are
c     something temporary, i.e., in the future some of this should come
c     in as input.
c
      ilimtyp = 1
      epsrden = 1.0d-10
c
c     If the user does not want any limiting, then don't compute
c     anything.
c
      if(ilimtyp.eq.0) then
         do 10 n=1,4
            phimins(n) = 1.0d0
            phiplus(n) = 1.0d0
   10    continue
         return
      end if
c
c     Compute the ratios of consecutive gradients (r_minus and r_plus).
c
      do 100 n=1,4
         aux1 = qpril1(n) - qpril2(n)
         aux2 = qprir2(n) - qprir1(n)
c
c        Check to avoid zero denominator for the "minus" side.
c
         if(dabs(aux1).ge.epsrden) then
c
c           OK, this is the general case, computer r_minus.
c
            rminus(n) = (qprir1(n) - qpril1(n)) / aux1
         else
c
c           Otherwise, I will assume that we are in a region of smooth
c           flow and turn the limiter off.
c
            rminus(n) = 1.0d0
         end if
c
c        Check to avoid zero denominator for the "plus" side.
c
         if(dabs(aux2).ge.epsrden) then
c
c           OK, this is the general case, computer r_plus.
c
            rplus(n) = (qprir1(n) - qpril1(n)) / aux2
         else
c
c           Otherwise, I will assume that we are in a region of smooth
c           flow and turn the limiter off.
c
            rplus(n) = 1.0d0
         end if
c
  100 continue
c
c     Computation of the "minmod" limiter.
c
      if(ilimtyp.eq.1) then
         do 200 n=1,4
c
c           Computation for the "minus" side.
c
            if(rminus(n).le.0.0d0) then
               phimins(n) = 0.0d0
            else
               if(rminus(n).ge.1.0d0) then
                  phimins(n) = 1.0d0
               else
                  phimins(n) = rminus(n)
               end if
            end if
c
c           Computation for the "plus" side.
c
            if(rplus(n).le.0.0d0) then
               phiplus(n) = 0.0d0
            else
               if(rplus(n).ge.1.0d0) then
                  phiplus(n) = 1.0d0
               else
                  phiplus(n) = rplus(n)
               end if
            end if
c
  200    continue
         return
      end if
c
c     If we get to this point, something is wrong. Print error message
c     and stop the code.
c
      write(*,901) ilimtyp
      stop
c
c     Formats.
c
  901 format(///,5x,'Subroutine LIMITER: You selected a limiter which',
     & ' has not been',/,5x,'implemented yet.',/,5x,'Valid options ',
     & 'are: (0) no limiter at all, and',/,24x,'(1) minmod.',/,5x,
     & 'You have selected option no. ',i6)
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine PRIMIV
c
c     This routine calculates the vector of primitive variables for
c     a specified volume (triangle).
c
c     The triangle (volume) is specified by its locator in the
c     "itable" list. Primitive variables are (p,u,v,T). However,
c     I prefer to use internal energy (e_i) instead of temperature (T).
c
c     The routine also checkes whether the specified volume is a ghost
c     volume. If it is, it will get the necessary information from the
c     appropriate array (in this case, "rghost" instead of "q").
c
c     Input variables are: itt     > identifies the desired triangle,
c                          qprimit > (p,u,v,e_i) for triangle itt.
c
      subroutine primiv(itt,qprimit)
c
      include 'Common_rm_v2.0.f'
c
      dimension qprimit(4)
c
c     First, we have to check if this is a ghost volume.
c
      if(itt.lt.0) go to 100
c
c     OK, we are dealing with a real volume. Just compute what you
c     need. Necessary information comes from the "q" array.
c
      aux1 = q(itt,4) - (0.5d0 * ((q(itt,2) * q(itt,2)) +
     &                            (q(itt,3) * q(itt,3))) / q(itt,1))
      qprimit(1) = (gamma - 1.0d0) * aux1
      qprimit(2) = q(itt,2) / q(itt,1)
      qprimit(3) = q(itt,3) / q(itt,1)
      qprimit(4) = aux1 / q(itt,1)
c
      return
c
c     We are dealing with a ghost volume in this case. Therefore,
c     information should come from the "rghost" array.
c
  100 continue
      ittp = - itt
      aux1 = rghost(ittp,4) - (0.5d0 * ((rghost(ittp,2)*rghost(ittp,2))
     &             + (rghost(ittp,3)*rghost(ittp,3))) / rghost(ittp,1))
      qprimit(1) = (gamma - 1.0d0) * aux1
      qprimit(2) = rghost(ittp,2) / rghost(ittp,1)
      qprimit(3) = rghost(ittp,3) / rghost(ittp,1)
      qprimit(4) = aux1 / rghost(ittp,1)
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine TESTINP
c
c     Given the coordinates of a point and the number of a triangle
c     (in the "itable" list), this routine returns whether the point
c     is within the triangle or not.
c     The flag is the variable "locflag" which has the following
c     values: locflag = 1 > point inside triangle;
c                       0 > point outside triangle.
c     Variable "locflag" always comes in with value = 1.
c     The other input variable are (obviously):
c             itt  >  identifies the triangle,
c             (xxx,yyy)  >  coordinates of the point.
c
      subroutine testinp(itt,locflag,xxx,yyy)
c
      include 'Common_rm_v2.0.f'
c
c     Testing performed here is based on the area of the triangles
c     formed by the new point and the edges of the existing triangle.
c     Since the nodes of the existing triangle are already ordered
c     in a counterclockwise fashion, all three areas will be positive
c     if point is inside triangle, and at least one of them will be
c     negative otherwise.
c
c     Get information on the nodes of the given triangle.
c
      inp1 = itable(itt,1)
      inp2 = itable(itt,2)
      inp3 = itable(itt,3)
c
c     Test triangle (inp1,inp2,new_point).
c
      xx1 = xynode(inp1,1)
      yy1 = xynode(inp1,2)
      xx2 = xynode(inp2,1)
      yy2 = xynode(inp2,2)
c
      call volumes(xx1,yy1,xx2,yy2,xxx,yyy,area)
c
      if(area.lt.0.0d0) then
         locflag = 0
         return
      end if
c
c     Test triangle (inp2,inp3,new_point).
c
      xx1 = xynode(inp2,1)
      yy1 = xynode(inp2,2)
      xx2 = xynode(inp3,1)
      yy2 = xynode(inp3,2)
c
      call volumes(xx1,yy1,xx2,yy2,xxx,yyy,area)
c
      if(area.lt.0.0d0) then
         locflag = 0
         return
      end if
c
c     Test triangle (inp3,inp1,new_point).
c
      xx1 = xynode(inp3,1)
      yy1 = xynode(inp3,2)
      xx2 = xynode(inp1,1)
      yy2 = xynode(inp1,2)
c
      call volumes(xx1,yy1,xx2,yy2,xxx,yyy,area)
c
      if(area.lt.0.0d0) locflag = 0
c
c     OK, the test is concluded. The calling routine will know what
c     to do.
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine VLFLUX2
c
c     This subroutine computes the E_plus, E_minus, F_plus and F_minus
c     flux vectors, according to the Van Leer flux vector splitting
c     scheme, for one interface each time it is called. A 2nd order
c     flux calculation is being used here, based on a MUSCL approach.
c
c     The procedure implemented performs linear reconstruction of
c     primitive variables, i.e., (p,u,v,T). Obviously, as usually,
c     I don't work with temperature (T), but with the internal
c     energy (e_i).
c
c     The routine actually already returns the complete interface flux
c     vectors to the calling routine. That is, it actually returns
c        E_ik = E_plus + E_minus     (in array "evec"), and
c        F_ik = F_plus + F_minus     (in array "fvec").
c     It turns out that these arrays ("evec" and "fvec") are the same
c     ones used by the central difference-type discretization (i.e.,
c     Jameson's method) in the present implementation. Don't get
c     confused !!! All I'm saying here is that the same Fortran arrays
c     are used in the two cases...
c
      subroutine vlflux2(j)
c
      include 'Common_rm_v2.0.f'
c
c     Some temporary storage.
c
      dimension qpril1(4), qpril2(4), qprir1(4), qprir2(4)
      dimension qpminus(4), qpplus(4)
      dimension phimins(4), phiplus(4)
c
c     Get the information from the edge-based database array.
c
      in1 = iedges(j,1)
      in2 = iedges(j,2)
      i = iedges(j,3)
      nb = iedges(j,4)
c
c     Decide who is left and right (for each direction).
c
      dxik = xynode(in2,1) - xynode(in1,1)
      dyik = xynode(in2,2) - xynode(in1,2)
c
      if(dyik.ge.0.0d0) then
         ixleft1 = i
         ixleft2 = iedaux(j,1)
         ixrigh1 = nb
         ixrigh2 = iedaux(j,2)
      else
         ixleft1 = nb
         ixleft2 = iedaux(j,2)
         ixrigh1 = i
         ixrigh2 = iedaux(j,1)
      end if
c
      if(dxik.le.0.0d0) then
         iyleft1 = i
         iyleft2 = iedaux(j,1)
         iyrigh1 = nb
         iyrigh2 = iedaux(j,2)
      else
         iyleft1 = nb
         iyleft2 = iedaux(j,2)
         iyrigh1 = i
         iyrigh2 = iedaux(j,1)
      end if
c
c     Calculation of the x-direction fluxes.
c
c     Make sure you have the right information.
c
      itleft1 = ixleft1
      itleft2 = ixleft2
      itrigh1 = ixrigh1
      itrigh2 = ixrigh2
c
c     Compute the necessary primitive variables.
c
      call primiv(itleft1,qpril1)
c
      if(itleft2.eq.itleft1) then
         do 60 idummy=1,4
            qpril2(idummy) = qpril1(idummy)
   60    continue
      else
         call primiv(itleft2,qpril2)
      end if
c
      call primiv(itrigh1,qprir1)
c
      if(itrigh2.eq.itrigh1) then
         do 70 idummy=1,4
            qprir2(idummy) = qprir1(idummy)
   70    continue
      else
         call primiv(itrigh2,qprir2)
      end if
c
c     Compute the limiter.
c
      call limiter(ilimtyp,qpril1,qpril2,qprir1,qprir2,phimins,phiplus)
c
c     Perform the linear reconstruction of left and right states.
c     (Remember, this is a "limited" linear reconstruction.)
c
      do 100 idummy=1,4
         qpminus(idummy) = qpril1(idummy) + 0.5d0 * phimins(idummy)
     &                   * (qpril1(idummy) - qpril2(idummy))
         qpplus(idummy) = qprir1(idummy) - 0.5d0 * phiplus(idummy)
     &                  * (qprir2(idummy) - qprir1(idummy))
  100 continue
c
c     Compute the fluxes in the x-direction.
c
      call vlfluxx(qpminus,qpplus)
c
c     Calculation of the y-direction fluxes.
c
c     First, check whether you really need to perform the
c     reconstruction (or if you could use the previous results).
c
      if(iyleft1.eq.itleft1 .and. iyrigh1.eq.itrigh1) go to 250
c
c     Otherwise, just make sure you have the right information and
c     repeat the procedure.
c
      itleft1 = iyleft1
      itleft2 = iyleft2
      itrigh1 = iyrigh1
      itrigh2 = iyrigh2
c
c     Compute the necessary primitive variables.
c
      call primiv(itleft1,qpril1)
c
      if(itleft2.eq.itleft1) then
         do 160 idummy=1,4
            qpril2(idummy) = qpril1(idummy)
  160    continue
      else
         call primiv(itleft2,qpril2)
      end if
c
      call primiv(itrigh1,qprir1)
c
      if(itrigh2.eq.itrigh1) then
         do 170 idummy=1,4
            qprir2(idummy) = qprir1(idummy)
  170    continue
      else
         call primiv(itrigh2,qprir2)
      end if
c
c     Compute the limiter.
c
      call limiter(ilimtyp,qpril1,qpril2,qprir1,qprir2,phimins,phiplus)
c
c     Perform the linear reconstruction of left and right states.
c     (Remember, this is a "limited" linear reconstruction.)
c
      do 200 idummy=1,4
         qpminus(idummy) = qpril1(idummy) + 0.5d0 * phimins(idummy)
     &                   * (qpril1(idummy) - qpril2(idummy))
         qpplus(idummy) = qprir1(idummy) - 0.5d0 * phiplus(idummy)
     &                  * (qprir2(idummy) - qprir1(idummy))
  200 continue
c
c     Compute the fluxes in the y-direction.
c
  250 continue
      call vlfluxy(qpminus,qpplus)
c
      return
c
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine VLFLUXX
c
c     This subroutine does the actual computation for obtaining the
c     E_plus and E_minus flux vectors for the Van Leer flux vector
c     splitting scheme. This routine assumes a MUSCL approach for the
c     definition of the flux vectors, and it also assumes that the
c     data is given in terms of primitive variables (p,u,v,ei).
c
c     A word of caution: Primitive variables are usually taken as
c     (p,u,v,T). I prefer, however, to work with the internal energy
c     (ei).
c
c     Another comment: This routine actually does not care whether
c     your computation is 1st or 2nd order (or even higher order).
c     All that it does is: given q_minus and q_plus, it computes the
c     corresponding E_plus(q_minus) and E_minus(q_plus), and sums
c     them up such that
c
c            E_ik = E_plus(q_minus) + E_minus(q_plus)
c
c
      subroutine vlfluxx(qpminus,qpplus)
c
      include 'Common_rm_v2.0.f'
c
      dimension eikplus(4), eikminus(4)
      dimension qpminus(4), qpplus(4)
c
c     Define some constants.
c
      const = gamma - 1.0d0
      const2 = (gamma * gamma) - 1.0d0
c
c     Perform the calculation for the "minus" side, i.e., for the
c     "left" side.
c
c     (Careful, however, because the meaning of left and right here
c     does not have anything to do with our usual convention of
c     orientation of an edge.)
c
c     Assign the usual names to the variables.
c
      p = qpminus(1)
      u = qpminus(2)
      v = qpminus(3)
c
c     Perform preliminary computations.
c
      rho = p / (const * qpminus(4))
      et = rho * (qpminus(4) + 0.5d0 * (u*u + v*v))
      a = dsqrt((gamma * p) / rho)
      amx = u / a
c
c     Compute the E_plus(q_minus) flux vector.
c
      if(amx.ge.1.0d0) then
c
c        For Mx > 1: E_plus = E.
c
         eikplus(1) = rho * u
         eikplus(2) = eikplus(1) * u + p
         eikplus(3) = eikplus(1) * v
         eikplus(4) = (et + p) * u
      else
c
c        For Mx < -1: E_plus = 0.
c
         if(amx.le.-1.0d0) then
            eikplus(1) = 0.0d0
            eikplus(2) = 0.0d0
            eikplus(3) = 0.0d0
            eikplus(4) = 0.0d0
c
         else
c
c        For -1 < Mx < 1: use the Van Leer formulas.
c
            fpmass = rho * a * (0.5d0 * (amx + 1.0d0))
     &                       * (0.5d0 * (amx + 1.0d0))
            aux1 = const * u + 2.0d0 * a
            eikplus(1) = fpmass
            eikplus(2) = fpmass * aux1 / gamma
            eikplus(3) = fpmass * v
            eikplus(4) = fpmass * (((aux1*aux1)/(2.0d0*const2))
     &                          + (0.5d0 * v * v))
         end if
      end if
c
c     Perform the calculation for the "plus" side, i.e., for the
c     "right" side.
c
c     (Again, the same comment applies here. Be careful because the
c     meaning of left and right here does not have anything to do with
c     our usual convention of orientation of an edge.)
c
c     Assign the usual names to the variables.
c
      p = qpplus(1)
      u = qpplus(2)
      v = qpplus(3)
c
c     Perform preliminary computations.
c
      rho = p / (const * qpplus(4))
      et = rho * (qpplus(4) + 0.5d0 * (u*u + v*v))
      a = dsqrt((gamma * p) / rho)
      amx = u / a
c
c     Compute the E_minus(q_plus) flux vector.
c
      if(amx.ge.1.0d0) then
c
c        For Mx > 1: E_minus = 0.
c
         eikminus(1) = 0.0d0
         eikminus(2) = 0.0d0
         eikminus(3) = 0.0d0
         eikminus(4) = 0.0d0
      else
c
c        For Mx < -1: E_minus = E.
c
         if(amx.le.-1.0d0) then
            eikminus(1) = rho * u
            eikminus(2) = eikminus(1) * u + p
            eikminus(3) = eikminus(1) * v
            eikminus(4) = (et + p) * u
c
         else
c
c        For -1 < Mx < 1: use the Van Leer formulas.
c
            fmmass = - rho * a * (0.5d0 * (amx - 1.0d0))
     &                         * (0.5d0 * (amx - 1.0d0))
            aux2 = const * u - 2.0d0 * a
            eikminus(1) = fmmass
            eikminus(2) = fmmass * aux2 / gamma
            eikminus(3) = fmmass * v
            eikminus(4) = fmmass * (((aux2*aux2)/(2.0d0*const2))
     &                           + (0.5d0 * v * v))
         end if
      end if
c
c     Finally, add the "plus" and "minus" fluxes in order to obtain
c     the flux vector at the interface (E_ik).
c
      evec(1) = eikplus(1) + eikminus(1)
      evec(2) = eikplus(2) + eikminus(2)
      evec(3) = eikplus(3) + eikminus(3)
      evec(4) = eikplus(4) + eikminus(4)
c
      return
      end
c
c ---------------------------------------------------------------------
c
c     Subroutine VLFLUXY
c
c     This subroutine does the actual computation for obtaining the
c     F_plus and F_minus flux vectors for the Van Leer flux vector
c     splitting scheme. This routine assumes a MUSCL approach for the
c     definition of the flux vectors, and it also assumes that the
c     data is given in terms of primitive variables (p,u,v,ei).
c
c     A word of caution: Primitive variables are usually taken as
c     (p,u,v,T). I prefer, however, to work with the internal energy
c     (ei).
c
c     Another comment: This routine actually does not care whether
c     your computation is 1st or 2nd order (or even higher order).
c     All that it does is: given q_minus and q_plus, it computes the
c     corresponding F_plus(q_minus) and F_minus(q_plus), and sums
c     them up such that
c
c            F_ik = F_plus(q_minus) + F_minus(q_plus)
c
c
      subroutine vlfluxy(qpminus,qpplus)
c
      include 'Common_rm_v2.0.f'
c
      dimension fikplus(4), fikminus(4)
      dimension qpminus(4), qpplus(4)
c
c     Define some constants.
c
      const = gamma - 1.0d0
      const2 = (gamma * gamma) - 1.0d0
c
c     Perform the calculation for the "minus" side, i.e., for the
c     "left" side.
c
c     (Careful, however, because the meaning of left and right here
c     does not have anything to do with our usual convention of
c     orientation of an edge.)
c
c     Assign the usual names to the variables.
c
      p = qpminus(1)
      u = qpminus(2)
      v = qpminus(3)
c
c     Perform preliminary computations.
c
      rho = p / (const * qpminus(4))
      et = rho * (qpminus(4) + 0.5d0 * (u*u + v*v))
      a = dsqrt((gamma * p) / rho)
      amy = v / a
c
c     Compute the F_plus(q_minus) flux vector.
c
      if(amy.ge.1.0d0) then
c
c        For My > 1: F_plus = F.
c
         fikplus(1) = rho * v
         fikplus(2) = fikplus(1) * u
         fikplus(3) = fikplus(1) * v + p
         fikplus(4) = (et + p) * v
      else
c
c        For My < -1: F_plus = 0.
c
         if(amy.le.-1.0d0) then
            fikplus(1) = 0.0d0
            fikplus(2) = 0.0d0
            fikplus(3) = 0.0d0
            fikplus(4) = 0.0d0
c
         else
c
c        For -1 < My < 1: use the Van Leer formulas.
c
            fpmass = rho * a * (0.5d0 * (amy + 1.0d0))
     &                       * (0.5d0 * (amy + 1.0d0))
            aux1 = const * v + 2.0d0 * a
            fikplus(1) = fpmass
            fikplus(2) = fpmass * u
            fikplus(3) = fpmass * aux1 / gamma
            fikplus(4) = fpmass * (((aux1*aux1)/(2.0d0*const2))
     &                          + (0.5d0 * u * u))
         end if
      end if
c
c     Perform the calculation for the "plus" side, i.e., for the
c     "right" side.
c
c     (Again, the same comment applies here. Be careful because the
c     meaning of left and right here does not have anything to do with
c     our usual convention of orientation of an edge.)
c
c     Assign the usual names to the variables.
c
      p = qpplus(1)
      u = qpplus(2)
      v = qpplus(3)
c
c     Perform preliminary computations.
c
      rho = p / (const * qpplus(4))
      et = rho * (qpplus(4) + 0.5d0 * (u*u + v*v))
      a = dsqrt((gamma * p) / rho)
      amy = v / a
c
c     Compute the F_minus(q_plus) flux vector.
c
      if(amy.ge.1.0d0) then
c
c        For My > 1: F_minus = 0.
c
         fikminus(1) = 0.0d0
         fikminus(2) = 0.0d0
         fikminus(3) = 0.0d0
         fikminus(4) = 0.0d0
      else
c
c        For My < -1: F_minus = F.
c
         if(amy.le.-1.0d0) then
            fikminus(1) = rho * v
            fikminus(2) = fikminus(1) * u
            fikminus(3) = fikminus(1) * v + p
            fikminus(4) = (et + p) * v
c
         else
c
c        For -1 < My < 1: use the Van Leer formulas.
c
            fmmass = - rho * a * (0.5d0 * (amy - 1.0d0))
     &                         * (0.5d0 * (amy - 1.0d0))
            aux2 = const * v - 2.0d0 * a
            fikminus(1) = fmmass
            fikminus(2) = fmmass * u
            fikminus(3) = fmmass * aux2 / gamma
            fikminus(4) = fmmass * (((aux2*aux2)/(2.0d0*const2))
     &                           + (0.5d0 * u * u))
         end if
      end if
c
c     Finally, add the "plus" and "minus" fluxes in order to obtain
c     the flux vector at the interface (F_ik).
c
      fvec(1) = fikplus(1) + fikminus(1)
      fvec(2) = fikplus(2) + fikminus(2)
      fvec(3) = fikplus(3) + fikminus(3)
      fvec(4) = fikplus(4) + fikminus(4)
c
      return
      end
c
c ---------------------------------------------------------------------
c
c       SUBROUTINE DynMesh                           Date: Mar./19/1993
c
c       This subroutine executes the movement of  an unstrucured mesh
c       following the algorithm suggested by Batina. Here we are using 
c       the file definitions in accord with the flow solver. These 
c       three files are necessary:
c
c          a) The file with the coordinates of the nodes (Fort.2);
c          b) The table with the identification of the nodes, orga-
c             nized by their position (Fort.12);
c          c) The table of neighborhood of the nodes (Fort.13).
c
c       These last two files are generated by the subroutine gen_f13
c       appended below.                                         
c       The following convention is used for the definition of the 
c       degree-of-freedom (idof):
c
c          (0) For the plunge displacement (h);
c          (1) For the pitch  displacement (alpha).
c
c       For the last displacement, we must give the position of the 
c       elastic axis (in percentage of chord, from the middle of the
c       profile, positive to the right).
c
c       And the following for the definition of the type of movement 
c       (imove):
c
c          (0) For harmonic oscillation;
c          (1) For exponentially shaped pulse;
c          (2) For step;
c          (3) For unit sample.                  
c
c       After we have displaced the mesh, we calculate the velocity of 
c       the nodes and the associated velocity of the grid cells.
c          
c     -----------------------------------------------------------------
c
      Subroutine dynmesh
c    
      include 'Common_rm_v2.0.f'
c
c     -----------------------------------------------------------------
c      Storing the coordinates into an intermediary array to calculate
c      the components of the grid velocity.
c     -----------------------------------------------------------------
c
      pi = dacos(-1.d0)
      if (niter.eq.1) then
         do 10 i=1,nnodes
            xhist(i,1) = x0(i)
            yhist(i,1) = y0(i)
            xhist(i,2) = x0(i)
            yhist(i,2) = y0(i)
   10    continue
      else         
         do 20 i=1,nnodes
            xhist(i,2) = xhist(i,1)
            yhist(i,2) = yhist(i,1)
            xhist(i,1) = xynode(i,1)
            yhist(i,1) = xynode(i,2)
   20    continue
      endif 
c
c     -----------------------------------------------------------------
c                   Moving the mesh.
c     -----------------------------------------------------------------
      time = niter * h                                                 
c
      if (imove.eq.0) then
        func = fsin(time)
      elseif (imove.eq.1) then
        func = fpulse(time)
      else
        func = 1.d0
      endif                              
c
      if (idof.eq.0) then
         if(imotion.eq.2 .and. imove.eq.1) then
            dht = ymax*func
            aux = ht
            ht = aux + dht
         else
            ht = ymax*func
         endif
         do 30 i=1,nnodes
c         
            xynode(i,1) = x0(i)
            xynode(i,2) = y0(i) + ht
c
   30    continue 
      else 
         if(imotion.eq.2 .and. imove.eq.1) then
            dalphat = alpham * func
            aux = alphat
            alphat = aux + dalphat
         else
            alphat = alpham * func
         endif
         do 40 i=1,nnodes 
c          
            xabs = x0(i) - xea
	      xynode(i,1) = xabs * dcos(alphat) + y0(i) * dsin(alphat) +
     &                    xea
	      xynode(i,2) = y0(i) * dcos(alphat) - xabs * dsin(alphat)
c
   40    continue
      endif
c
      return
c
      end
c
c     -----------------------------------------------------------------
c     Subroutine Velocity                           Noll(Feb.-06-2006)
c
c     This subroutine was part of the dynmesh subroutine. I separated it
c     because I wanted to call only this part in other point of the code.
c     -----------------------------------------------------------------
      Subroutine velocity
c    
      include 'Common_rm_v2.0.f'
      dimension vnode(nmaxnd,2)
c
c     -----------------------------------------------------------------
c            Calculating the components of the grid velocity.             
c     -----------------------------------------------------------------
      if (ivmesh.lt.2) then
         do 10 i=1,nnodes
            vnode(i,1) = (xynode(i,1) - xhist(i,1))/h
            vnode(i,2) = (xynode(i,2) - yhist(i,1))/h
   10    continue
      else
         do 20 i=1,nnodes
            vnode(i,1) = 0.5d0/h*(3.0d0*xynode(i,1) - 4.0d0*xhist(i,1)
     &                                              + xhist(i,2))
            vnode(i,2) = 0.5d0/h*(3.0d0*xynode(i,2) - 4.0d0*yhist(i,1)
     &                                            + yhist(i,2))
   20    continue
      endif
c
c     ----------------------------------------------------------------- 
c              Computing the componentes of the edge velocity.
c     -----------------------------------------------------------------
c                                               
      do 30 i=1,nsides
         in1 = iedges(i,1)
         in2 = iedges(i,2)
         nb = iedges(i,4)
         vside(i,1) = 0.5d0 * (vnode(in1,1) + vnode(in2,1))
         vside(i,2) = 0.5d0 * (vnode(in1,2) + vnode(in2,2))
         if (nb.lt.0) then
            uvolg(-nb) = vside(i,1)
            vvolg(-nb) = vside(i,2)
         endif
   30 continue    
c
      return
      end         
c
c     -----------------------------------------------------------------
c
c      This function calculates the value of the harmonic function.
c
c     -----------------------------------------------------------------
c
      Real*8 Function fsin(time)
c
      include 'Common_rm_v2.0.f'
c
      fsin = sin(akc*fsmach*time)     
c
c     -----------------------------------------------------------------
      return
      end                                         
c                                                                      
c     -----------------------------------------------------------------
c
c      This function calculates the value of the exponentially shaped 
c      pulse.
c
c     -----------------------------------------------------------------
c
      Real*8 Function fpulse(time)
c
      include 'Common_rm_v2.0.f'
c
      tb = time/tmax
c
      if (tb.lt.1.0d0) then
        cte1 = 2.d0 - 1.d0/(1.d0-tb)
        fpulse= 4.d0 * tb*tb * exp(cte1)
      else
        fpulse = 0.d0
      endif
c
c     -----------------------------------------------------------------
      return
      end                                         
c
c     *****************************************************************
c
c       This subroutine generates the conectivity table of nodes used
c       by the dynamic mesh algorithm for unstructured meshes.
c       Here only interior nodes are considered.
c       Three files are necessary:
c
c             a) The conectivity table (Fort.3);
c             b) The neighborhood table (Fort.4);
c             c) The table of nodes organized by their position
c                (Fort.12).
c
c
c     *****************************************************************
c
      Subroutine gen_f13
c
      include 'Common_rm_v2.0.f'
c
      open(12,file='Fort.12')
      open(13,file='Fort.13')
c     ----------------------------------------------------------------- 
c     Generating the intermediary array..
c     ----------------------------------------------------------------
      nwn=0
      nobn=0
      nin=0
c
      do 30 i=1,nnodes
        iflag=0
c
c       Boundary nodes.
c       ---------------
        do 10 jj=1,nvolg
          if (iflag.eq.0) then
c
            if (ighost(jj,3).eq.i .or. ighost(jj,4).eq.i) then
c
              if (ighost(jj,3).eq.i) then
                itab=ighost(jj,3)
              else
                itab=ighost(jj,4)
              endif
c
              if (ighost(jj,1).eq.1) then
                nwn=nwn+1
                iwn(nwn)=itab
              else
                nobn=nobn+1
                iobn(nobn)=itab
              endif
              iflag=1
c
            endif
c
          endif
   10   continue
c
c       Interior nodes.
c       ---------------
        do 20 jj=1,nvolg
          if (iflag.eq.0) then
c
            if (ighost(jj,3).ne.i .and. ighost(jj,4).ne.i) then
              nin=nin+1
              iin(nin)=i
              iflag=1
            endif
c
          endif
   20   continue
c
   30 continue
c
c     -----------------------------------------------------------------
c     Now we begin the loop to cover all the interior nodes.
c     -----------------------------------------------------------------
      do 60 ii=1,nin                      
c
        inode = iin(ii)
c
c       Searching for the first triangle that contains this node.
c       --------------------------------------------------------
        iflag = 0
        n = 1        
   40   if (iflag.eq.0) then
c       
          it = n
          in1=itable(n,1)
          in2=itable(n,2)
          in3=itable(n,3)
c
          if ((in1.eq.inode).or.(in2.eq.inode).or.(in3.eq.inode)) then
            itinit = it
            iflag = 1
          endif
          n = n+1
          goto 40
c
        endif
c
c       Searching for the first neighbor triangle and the first neighbor
c       node.
c       ----------------------------------------------------------------
        icount = 1
        call Neighb
c
c       Now we begin to look for the other neighbors until we come back
c       to the first one.
c       ---------------------------------------------------------------
   50   if (itneigh.ne.itinit) then
c
          it = itneigh
          icount = icount+1
c                                                                      
          in1=itable(it,1)
          in2=itable(it,2)
          in3=itable(it,3)
c
          call Neighb
c
          goto 50
c
        endif       
                   
c
        jtot(inode) = icount                    
        write(13,*) inode,icount,(neigh(inode,k),k=1,icount)
c                
   60 continue

c
c     -----------------------------------------------------------------
c     Writing into the file.
c     -----------------------------------------------------------------
c
      write(12,*) nwn,nobn,nin
      do 70 i=1,nwn
        write(12,*) iwn(i)
   70 continue
c
      do 80 i=1,nobn
        write(12,*) iobn(i)
   80 continue
c
      do 90 i=1,nin
        write(12,*) iin(i)
   90 continue 
c
      close(12)
      close(13)    
c     -----------------------------------------------------------------
      return
      end
c
c     -----------------------------------------------------------------
c
c     Being given a triangle that contains a determined node, this sub-
c     routine searchs for its neighbor and finds the farmost right nei-
c     ghbor of this particular node.
c
c     -----------------------------------------------------------------
c
      Subroutine Neighb
c
      include 'Common_rm_v2.0.f'
c     -----------------------------------------------------------------
c     Determining the neighbors.
c                  
      nb1=neighbor(it,1)
      nb2=neighbor(it,2)
      nb3=neighbor(it,3)
c
      if (inode.eq.in1) then
        itneigh = nb3
        i = in3
      endif
c
      if (inode.eq.in2) then
        itneigh = nb1
        i = in1
      endif
c
      if (inode.eq.in3) then
        itneigh = nb2
        i = in2
      endif              
      neigh(inode,icount) = i
c
c     -----------------------------------------------------------------
      return
      end
c     -----------------------------------------------------------------
c
c     This subroutine calculates the pressure distribution along the
c     wall boundaries and output the results in a convinient nondimen-
c     sionalized manner. The output is also order in a fashion it eases
c     the plotting using Tecplot.
c
c     -----------------------------------------------------------------
c
      Subroutine pressure    
c     
      include 'Common_rm_v2.0.f'
c
c     Define some temporary storage.
c
      common p(2,nmaxgh),x(2,nmaxgh),y(2,nmaxgh),
     &       dens(2,nmaxgh),entr(2,nmaxgh)
      common iup,ilow
c
c     Identify "ghost" volumes which are actually at a wall and wether they
c     are an upper or lower boundary. Also compute the pressure.
c
      iup = 0
      ilow = 0
      do 50 ig=1,nvolg
         ibtype = ighost(ig,1)
         if(ibtype.eq.1) then
            const = gamma - 1.0d0
            pwall =  const * (rghost(ig,4) - 0.5d0*
     &                    (rghost(ig,2)*rghost(ig,2) +
     &                     rghost(ig,3)*rghost(ig,3))/
     &                     rghost(ig,1))
            in1 = ighost(ig,3)
            in2 = ighost(ig,4)
            x1 = xynode(in1,1)
	    y1 = xynode(in1,2)
            x2 = xynode(in2,1)
            y2 = xynode(in2,2)
            if(x2.gt.x1) then
	       if(iflow.eq.1) then
                  iflag = 0
	       else
	          iflag = 1
	       end if
	    else
	       if(iflow.eq.1) then
                  iflag = 1
	       else
	          iflag = 0
	       end if
	    endif
	    if(iflag.eq.0) then
	       iup = iup + 1
	       p(1,iup) = pwall
               dens(1,iup) = rghost(ig,1)
	       x(1,iup) = (x1+x2)/2
	       y(1,iup) = (y1+y2)/2
	    else
               ilow = ilow + 1
	       p(2,ilow) = pwall
               dens(2,ilow) = rghost(ig,1)
	       x(2,ilow) = (x1+x2)/2
	       y(2,ilow) = (y1+y2)/2
	    endif
         endif
   50 continue
c  
c     Order the wall points
c
      do 70 i=1,iup
         do 80 ii=iup,i+1,-1
            if(x(1,i).gt.x(1,ii)) then
	       exchg = x(1,i)
	       x(1,i) = x(1,ii)
               x(1,ii) = exchg
	       exchg = y(1,i)
	       y(1,i) = y(1,ii)
               y(1,ii) = exchg
	       exchg = p(1,i)
	       p(1,i) = p(1,ii)
               p(1,ii) = exchg
	       exchg = dens(1,i)
	       dens(1,i) = dens(1,ii)
               dens(1,ii) = exchg
	    endif
   80    continue
   70 continue
c
      do 90 i=1,ilow
         do 100 ii=ilow,i+1,-1
            if(x(2,i).lt.x(2,ii)) then
	       exchg = x(2,i)
	       x(2,i) = x(2,ii)
               x(2,ii) = exchg
	       exchg = y(2,i)
	       y(2,i) = y(2,ii)
               y(2,ii) = exchg
	       exchg = p(2,i)
	       p(2,i) = p(2,ii)
               p(2,ii) = exchg
	       exchg = dens(2,i)
	       dens(2,i) = dens(2,ii)
               dens(2,ii) = exchg
	    endif
  100    continue
   90 continue
c
c     Compute pressure distribution at the wall. For external flow
c     cases, the quantity computed is the pressure coefficient Cp.
c     For internal flow cases with supersonic entrance, the quantity
c     computed is the ratio p/p_infty. Finally, for internal flow
c     cases with a subsonic entrance, the quantity computed is the
c     ratio p/P_t.
c
      if(iflow.eq.1) then
         const = 2.0d0 / (gamma * fsmach * fsmach)
         do 130 i=1,iup
            entr(1,i) = gamma*p(1,i) / (dens(1,i)**gamma) - 1.0d0
            p(1,i) = const * (gamma*p(1,i) - 1.0d0)
  130    continue
         do 135 i=1,ilow
            entr(2,i) = gamma*p(2,i) / (dens(2,i)**gamma) - 1.0d0
            p(2,i) = const * (gamma*p(2,i) - 1.0d0)
  135    continue
      elseif(iflow.eq.2) then
         do 140 i=1,iup
            p(1,i) = gamma * p(1,i)
  140    continue
         do 145 i=1,ilow
            p(2,i) = gamma * p(2,i)
  145    continue
      elseif(iflow.eq.3) then
         const = (2.0d0*gamma) / (gamma + 1.0d0)
         do 150 i=1,iup
            p(1,i) = const * p(1,i)
  150    continue
         do 155 i=1,ilow
            p(2,i) = const * p(2,i)
  155    continue
      endif
c          
c     Print results on the appropriate file
c     (steady,unit = 24; unesteady units=19).
c
      if(iflow2.eq.1) then 
        pi = acos(-1.0d0)
        alphai = alphat*180.d0/pi
        write(19,923)
	write(19,924) niter,alphai
c
        do 160 i=1,iup 
           x(1,i) = x(1,i) / cos(alphat)
           if (imove.eq.0) then
              y(1,i) = y(1,i) + x(1,i)*sin(alphat)
           else 
              y(1,i) = y(1,i) - ht
           endif
           write(19,905) x(1,i),y(1,i),-p(1,i)
  160   continue
c
        do 165 i=1,ilow
           x(2,i) = x(2,i) / cos(alphat)
           if (imove.eq.0) then
              y(2,i) = y(2,i) + x(2,i)*sin(alphat)
           else 
              y(2,i) = y(2,i) - ht
           endif
           write(19,905) x(2,i),y(2,i),-p(2,i)
  165   continue
c
        return
      endif
c
      if(iflow.eq.1) then
         write(20,901)
      elseif(iflow.eq.2) then
         write(20,902)
      else      
         write(20,903)
      end if
      write(20,904)
      if(iflow.eq.1) then
         do 200 i=1,iup
	    write(20,906) x(1,i),y(1,i),-p(1,i),entr(1,i)
  200    continue
         write(20,904)
	 do 210 i=1,ilow
            write(20,906) x(2,i),y(2,i),-p(2,i),entr(2,i)
  210    continue
      else
         do 220 i=1,iup
            write(20,905) x(1,i),y(1,i),p(1,i)
  220    continue
         write(20,904)
	 do 230 i=1,ilow
	    write(20,905) x(2,i),y(2,i),p(2,i)
  230    continue
      end if
c
c     -----------------------------------------------------------------
      return
c     -----------------------------------------------------------------
c
c     Formats
c
  901 format('VARIABLES = "X","Y","-Cp","S"')
  902 format('VARIABLES = "X","Y","P/Pinf"')
  903 format('VARIABLES = "X","Y","P/Pt"')
  904 format(/,'ZONE')
  905 format(3(f15.8,5x))
  906 format(4(f15.8,5x))
  923 format('VARIABLES = "X","Y","-Cp"')
  924 format('ZONE T="It = ',i6,1x,'Alpha = ',e9.3,'"') 

      end
c     -----------------------------------------------------------------
