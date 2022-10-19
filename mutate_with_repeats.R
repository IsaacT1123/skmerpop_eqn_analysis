simulation <- function(r=c(1,0.0,0.0), L=100000, d=0.01, rho=0.0) {
  print(data.frame(r1=r[1],r2=r[2],r3=r[3],L=L,d=d,rho=rho));
  p=1-(1-d)^31
  
  G1=c(rep(seq(1,10000,length=L*r[1]),1),
       rep(seq(20000,29999,length=L*r[2]/2),2),
       rep(seq(30000,39999,length=L*r[3]/3),3))
  #print(length(G1))
  
  tryCatch( {
    #create distribution for mutation selection
    mutcount = rpois(n = 1,lambda=L*p)
    #select positions to not be mutated
    nonmut_pos = sample(1:L,L-mutcount)
    
    #select positions which will be mutate into existing kmers
    mut_pos_silent = sample(1:L,  rho*mutcount )
    
    #select positions that will mutate into unique kmers
    mut_pos_new = seq(90000,99999,length= (1-rho)*mutcount )
    
    G2= c(G1[nonmut_pos],G1[mut_pos_silent],mut_pos_new)
    #print(c(length(G2),length(nonmut_pos),length(mut_pos_new),length(mut_pos_silent)))
    
    G2r = table(table(G2))
    G1r = table(table(G1))
    #print(list(r*L,G1r, G2r))
    
    ed = function(G1,G2) {
      I = intersect(G1,G2); 
      U=union(G1,G2); 
      j=length(I)/length(U);
      data.frame(
        rho=rho,
        d=d,
        r1= r[1],
        r2= r[2],
        r3 = r[3],
        trued=d,
        pm = p,
        sumr1 = sum(G1r),
        sumr2 = sum(G2r),
        I_real=length(I),
        I_skmer=L*(1-p),
        I_raR1=sum(sapply(1:length(G1r), function(i) G1r[i]* ( 1 - ((1-rho)*p)^i))),
        I_raR2=sum(sapply(1:length(G2r), function(i) G2r[i]* ( 1 - ((1-rho)*p)^i))),
        d_skmerest=1-(2*j/(1+j))^(1/31));}
    
    ed(G1,G2)}, error= function(error) {
      data.frame( rho=rho,
                  d=d,
                  r1= r[1],
                  r2= r[2],
                  r3 = r[3],
                  trued=d,
                  pm = p,
                  sumr1 = NA,
                  sumr2 = NA,
                  I_real=NA,
                  I_skmer=NA,
                  I_raR1=NA,
                  I_raR2=NA,
                  d_skmerest=NA)
    })
}
​
sim <- function(r2, L, d, rho, rep=2) simulation(r=c(1-3/2*r2,r2,r2/2),L=L, d=d, rho=rho)

res=do.call("rbind", 
            apply(expand.grid(r2=c(0,0.2,0.4), L=100000, d=c(0.01,0.03), rho=c(0.0,0.2,1)), 1, 
                  function(x) do.call(sim,as.list(x))) )
res