# This code needs more testing. However, it looks like the 'dry' dam-break
# does conserve energy, while the 'wet' dam break does not. Theoretically I 
# think this is due to the shock in the latter case.

# In Steve's notes, the chapter on the dam-break solution mentions a solution
# for the Shock height, and how this is not perfectly related to the energy
# dissipation -- and refers to a missing figure!
#

#
# Does the dam-break solution conserve energy? This dry dam-break seems to
#
dry_dam_break<-function(t=2, L=100, H1=1, N=1001, G=9.8){
    
    if(L < sqrt(G*H1)*t) stop('L is not long enough to support this t')

    c0 = sqrt(G * H1)
    x = seq(-L, L, len=N)
    
    U = x*0
    H = x*0
    
    R1 = which((x > -2*c0*t) & (x < c0*t))
    H[R1] = 1/(9*G)*(x[R1]/t + 2*c0)^2
    U[R1] = 2/3*(x[R1]/t - c0)
    
    H[(x > c0*t)] = H1

    energy = (x[2]-x[1])*sum( U^2*H + G*H^2) * 0.5
    #print(energy)
    
    out = list(x=x, h=H, u=U, energy=energy)

    return(out)

}

# Run a few cases -- yes energy is conserved
#tmp = dry_dam_break(0.0001, N=10001)
#tmp = dry_dam_break(1, N=10001)
#tmp = dry_dam_break(5, N=10001)

#
# What about the wet dam break? This needs more testing. The energy plots
# made later on clearly have some artefacts
#
wet_dam_break<-function(t=2, L=100, H0=0.1, H1=1, N=1001, G=9.8){

    if(H0 > H1) stop('H0 must be < H1')

    if(L < sqrt(G*H1)*t) stop('L is not long enough to support this t')

    #To calculate the Stoker solution, we minimize this function, Wu et al., (1999)
    allvars<-function(v){
        #u2 is velocity in shock zone
        u2=v[1]
        #h2 is depth
        h2=v[2]
        #epdot is wave speed
        epdot=v[3]

        if(h2 < H0) return(9e+10)

        #Preliminary celerities
        C0=sqrt(G*H0)
        C1=sqrt(G*H1)
        C2=sqrt(G*h2)

        # See eqns 44 - 47 in Wu et al. (1999) - these should all be zero for the solution
        eq1= epdot/C0 - 1/4*C0/epdot*(1+sqrt(1+8*(epdot/C0)^2)) -u2/C0
        eq2= 1/sqrt(2)*(sqrt(1+8*(epdot/C0)^2)-1  )^0.5 -C2/C0
        eq3= u2+2*(C2-C1)

        return(abs(eq1)+abs(eq2)+abs(eq3))

    }

    # Zoppou and Roberts solution
    # In my PhD code, I suggested there were problems with this in some cases.
    # So the values were used as initial guesses for another optimization
    # This seems to work
    S2fun<-function(S2){
        abs(-S2+2*sqrt(G*H1)+G/(4*S2)*(H0+sqrt(H0^2+8*H0*S2^2/(G)) ) -(2*G*sqrt(H0^2+8*H0*S2^2/(G)) -2*G*H0 )^0.5)
    }
    a11 = optimize(S2fun, interval=c(1.0E-9,1000),tol=1E-12)
    S2 = a11$minimum
    # So now we have the shock speed. We need u2, and h2
    u2 = S2-G/(4*S2)*(H0+sqrt(H0^2+8*H0*S2^2/(G)))

    # h2 also needs minimizing --- however, I think this gives incorrect values
    h2fun<-function(h2){
        #abs(-h2 + h0/2*sqrt(1+8*(2*h2/(h2-h0)*(sqrt(h1)-sqrt(h2))/sqrt(h0) )^2  ) -1/2)
        abs(-h2 + 1/2*sqrt(H0^2+8*H0*(2*h2/(h2-H0)*(sqrt(H1)-sqrt(h2)) )^2  ) -1/2)
    }
    a11 = optimize(h2fun, interval=c(H0+1E-9,H1)) #This interval is important, should be sensible I guess
    h2 = a11$minimum

    # Find the minimum of allvars- practically, this should set all equations to 0.
    # For initial values, make these guesses
    #S2 = sqrt(G*H0)
    #h2 = 0.5*(H1 + H0)
    #u2 = sqrt(G*h2)

    #a11 = optim(c(u2,h2,S2),allvars,lower=c(0,H0,0),method='L-BFGS-B', control=list(abstol=1e-12))
    a11 = optim(c(u2,h2,S2), allvars, control=list(abstol=1e-15, reltol=1e-12))
    u2 = a11$par[1]
    h2 = a11$par[2]
    S2 = a11$par[3]

    # x values
    xz = seq(-L,L, len=N)

    # Pack the solution
    u = rep(0, N)
    h = rep(H1,N)
    b = which( (xz> -t*sqrt(G*H1))&(xz <= t*(u2-sqrt(G*h2)) ) )
    u[b] = 2/3*(sqrt(G*H1)+xz[b]/t)
    h[b] = 4/(9*G)*(sqrt(G*H1)-xz[b]/(2*t) )^2
    #Now find the next interval
    b = which( ( xz>t*(u2-sqrt(G*h2)))&(xz < t*S2 ))
    u[b] = u2
    h[b] = h2
    #And the final interval
    b = which(xz>=t*S2)
    u[b] = 0
    h[b] = H0

    energy = (xz[2]-xz[1])*sum( u^2*h + G*h^2) * 0.5
    #print(energy)

    out = list(x=xz, h=h, u=u, energy=energy)
    return(out)
}

# Here we see the wet dam break solution is bleeding energy over time,
# (assuming I correctly coded it!)
#tmp = wet_dam_break(0.0001, N=10001)
#tmp = wet_dam_break(1, N=10001)
#tmp = wet_dam_break(5, N=10001)
#tmp = wet_dam_break(10, N=10001)
#tmp = wet_dam_break(20, N=10001)

# Figure showing energy loss issues.
make_energy_loss_figure<-function(){
    png('Dam_break_energy_loss.png', width=12, height=8, units='in', res=300)
    par(mfrow=c(2,2))

    # Look at energy decay over time. This seems convergent in N
    ts = seq(0,15)
    energy_over_time = sapply(ts, f<-function(t) wet_dam_break(t, N=100001)$energy)
    plot(ts, energy_over_time, main='Energy decay of wet dam-break solution', 
         xlab='Time', ylab='Energy')

    # What about as we get closer to a dry dam-break?
    energy_over_time_B = sapply(ts, f<-function(t) wet_dam_break(t, N=100001, H0=1e-10)$energy)
    plot(ts, energy_over_time_B, main='Energy in the limit of a dry dam-break solution', 
         xlab='Time', ylab='Energy')

    sol1 = wet_dam_break(t=15, H0=1e-10)
    sol2 = dry_dam_break(t=15)
    plot(sol1$x, sol1$h, t='l', 
         main='Test: dry dam-break vs wet-dam-break with H0 --> 0')
    points(-sol2$x, sol2$h, t='l', col='red')


    # Get energy loss over time as a function of H0
    H0s = seq(0.01, 0.99, by=0.01)
    energy_loss_rate_per_sec<-function(H0){

        # Compute over 10 seconds
        e0 = wet_dam_break(0, N=10001, H0=H0)$energy
        e1 = wet_dam_break(10, N=10001, H0=H0)$energy

        return((e1-e0)/(10))
    }

    energy_loss = sapply(H0s, energy_loss_rate_per_sec)
    plot(H0s, energy_loss, main='Wet dam-break energy loss vs downstream depth \n FIXME Some convergence issues?', 
         xlab='Depth downstream', ylab='Energy loss / second')
    dev.off()
}
