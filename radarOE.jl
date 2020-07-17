using PyCall

pickle=pyimport("pickle")

f=pybuiltin("open")("dataVN_3D.pklz","rb")
dataOut=pickle["load"](f)

clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,
pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,
zKuL,zKaL,cmbSfcRainL,gvSfcRainL,cmbPWC=dataOut;

include("scatTables.jl")
using .scatTables

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
extW_r,scatW_r,gW_r,vfallW_r=scatTables.init()

xr=pyimport("xarray")
function saveTables()
    nwTX=xr["DataArray"](nwT)
    pwcTX=xr["DataArray"](pwcT)
    attKaTX=xr["DataArray"](attKaT)
    attKuTX=xr["DataArray"](attKuT)
    zKaTX=xr["DataArray"](zKaT)
    zKuTX=xr["DataArray"](zKuT)
    zKaTsX=xr["DataArray"](zKaTs)
    zKuTsX=xr["DataArray"](zKuTs)
    dmTsX=xr["DataArray"](dmTs)
    dmTX=xr["DataArray"](dmT)
    pwcTsX=xr["DataArray"](pwcTs)
    attKaTsX=xr["DataArray"](attKaTs)
    attKuTsX=xr["DataArray"](attKuTs)
    attWTsX=xr["DataArray"](attWTs)
    zWTsX=xr["DataArray"](zWTs)
    nwTsX=xr["DataArray"](nwTs)
    rateX=xr["DataArray"](rateT)
    ratesX=xr["DataArray"](rateTs)
    d=Dict("nwT"=>nwTX,"pwcT"=>pwcTX,"attKaT"=>attKaTX,
           "attKuT"=>attKuTX, "dmT"=>dmTX, "zKuT"=>zKuTX, "zKaT"=>zKaTX,
           "nwTs"=>nwTsX,"pwcTs"=>pwcTsX,"attKaTs"=>attKaTsX,
           "attKuTs"=>attKuTsX, "dmTs"=>dmTsX, "zKuTs"=>zKuTsX, "zKaTs"=>zKaTsX,
           "zWTs"=>zWTsX,"attWTs"=>attWTsX, "rate"=>rateX,
           "rateS"=>ratesX)
    tables=xr["Dataset"](d)
    tables["to_netcdf"]("tables_nsv_rho04_dmTs11.nc")
end

iread=0
if iread==1
    ns=9
    scatTables.getDmNwSF(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,ns,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    scatTables.getDmNwR(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    saveTables()
    exit(1)
end  
n4=pyimport("netCDF4")
fh=n4["Dataset"]("tables_nsv_rho06_dmTs11.nc")
zKuT=fh["variables"]["get"]("zKuT")["_get"]([0],[240],[1])
zKaT=fh["variables"]["get"]("zKaT")["_get"]([0],[240],[1])
attKuT=fh["variables"]["get"]("attKuT")["_get"]([0],[240],[1])
attKaT=fh["variables"]["get"]("attKaT")["_get"]([0],[240],[1])
dmT=fh["variables"]["get"]("dmT")["_get"]([0],[240],[1])
nwT=fh["variables"]["get"]("nwT")["_get"]([0],[240],[1])
pwcT=fh["variables"]["get"]("pwcT")["_get"]([0],[240],[1])
zKuTs=fh["variables"]["get"]("zKuTs")["_get"]([0],[240],[1])
zKaTs=fh["variables"]["get"]("zKaTs")["_get"]([0],[240],[1])
attKuTs=fh["variables"]["get"]("attKuTs")["_get"]([0],[240],[1])
attKaTs=fh["variables"]["get"]("attKaTs")["_get"]([0],[240],[1])
dmTs=fh["variables"]["get"]("dmTs")["_get"]([0],[240],[1])
nwTs=fh["variables"]["get"]("nwTs")["_get"]([0],[240],[1])
pwcTs=fh["variables"]["get"]("pwcTs")["_get"]([0],[240],[1])
attWTs=fh["variables"]["get"]("attWTs")["_get"]([0],[240],[1])
zWTs=fh["variables"]["get"]("zWTs")["_get"]([0],[240],[1])
rateTs=fh["variables"]["get"]("rateS")["_get"]([0],[240],[1])
rateT=fh["variables"]["get"]("rate")["_get"]([0],[240],[1])
include("attCorrection.jl")

push!(PyVector(pyimport("sys")["path"]), "./")

np=pyimport("numpy")
h=0.125*((176:-1:100).-100);
include("radarMod.jl")
fractS=Spline1D([0.,1.25,3.85,6.0,6.35,10.],[1,1,1,0.1,0,0],k=1)
fractD=Spline1D([0.,1.25,3.85,6.35,10.],[1.8,1.8,1.7,1.6,0.7])
h1=(176:-1:1)*0.125
fract=fractS(h1[100:165])

using PyPlot
using LinearAlgebra

fP=pybuiltin("open")("dmCov.pklz","rb")
#
s=pickle["load"](fP)
eVal,eVect,dmM=s
print(length(s))


h1node=h1[[100,115,125,145,155,165]]
scipy=pyimport("scipy.ndimage") 

#include("radarMod.jl")
using .radarMod
radarMod.setRadar_dr(0.125,66)
radarMod.setRainFract(fract)
radarMod.setdNw(zeros(66))
function simZ(Dm)
    drk=Δr
    fract=rainFract
    nz=size(dmInt)[1]
    piaKu=0
    zSim=zeros(nz)
    rrateL=zeros(nz)
    noise=randn(nz)*1.0
    noise=scipy["gaussian_filter"](noise,5)
    for i=1:nz
        if fract[i]<0.5
            Dm[i]*=0.8
        else
            Dm[i]*=1-0.4*(1-fract[i])
        end
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        att=att*10^noise[i]
        piaKu=piaKu+att*drk
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        piaKu+=att*drk
    end
    return piaKu,rrateL,zSim
end

function simZ(Dm)
    drk=Δr
    fract=rainFract
    nz=size(dmInt)[1]
    piaKu=0
    zSim=zeros(nz)
    rrateL=zeros(nz)
    noise=dNw
    for i=1:nz
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        att=att*10^noise[i]
        piaKu=piaKu+att*drk
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        piaKu+=att*drk
    end
    return piaKu,rrateL,zSim
end

function simZ_g(Dm)
    drk=Δr
    fract=rainFract
    nz=size(dmInt)[1]
    piaKu=0
    zSim=zeros(nz)
    rrateL=zeros(nz)
    noise=dNw
    dZdDm=zeros(nz,nz)
    for i=1:nz
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        if n1==240
            n1=239
        end
        if n1s==240
            n1s=239
        end
       
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        n2s1=n2s
        n21=n2
        if n2s==1
            n2s1=2
        end
        if n2==1
            n21=1
        end
         dDm=(1-fract[i])*dmTs[n2s1]+fract[i]*dmT[n21]-
        ((1-fract[i])*dmTs[n1s]+fract[i]*dmT[n1])+1e-3
        att2=(1-fract[i])*attKuTs[n2s1]+fract[i]*attKuT[n21]
        att=att*10^noise[i]
        att2=att2*10^noise[i]
        piaKu=piaKu+att*drk
        dpiaKu1=(att2-att)/dDm*Δr
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        zSim2=log10((1-fract[i])*10^(0.1*zKuTs[n2s1])+(fract[i])*10^(0.1*zKuT[n21]))*10.0-piaKu+10*noise[i]
        dZdDm1=(zSim2-zSim[i])/dDm
        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        piaKu+=att*drk
        dZdDm[i,i]=dZdDm[i,i]+dZdDm1-dpiaKu1
        if i+1<=nz
            dZdDm[i,i+1:end]=dZdDm[i,i+1:end].-2*dpiaKu1
        end
    end
    return piaKu,rrateL,zSim,dZdDm
end


eVal1=eVal.*2.0




nT=3000
zSimL=zeros(3*nT,66)
pRateL=zeros(3*nT,66)

eVal1=eVal.*1.0

dmL=zeros(nT,66)
for i=1:nT
    global dmP=transpose(randn(1)*eVal1[6].*eVect[6,:]+randn(1)*eVal1[5].*eVect[5,:]+randn(1)*eVal1[4].*eVect[4,:]+
    randn(1)*eVal1[3].*eVect[3,:]+randn(1)*eVal1[2].*eVect[2,:]+randn(1)*eVal1[1].*eVect[1,:])+dmM
    dmS=Spline1D(h1node[4:-1:1],dmP[4:-1:1],k=2)
    global dmInt=dmS(h1[100:165]).+randn(66)*0.2
    nz=66
    noise=randn(nz)*1.0
    noise=scipy["gaussian_filter"](noise,10)
    dmInt=dmInt.+0.5*noise
    for k=1:66
        if fract[k]<0.5
            dmInt[k]*=0.8
        else
            dmInt[k]*=1-0.4*(1-fract[k])
        end
    end
    piaKu,rrate,zSim=simZ(dmInt)
    zSimL[i,:]=zSim
    #pRateL[i,:]=rrate
    dmL[i,:]=dmInt
end

using Statistics
dmMV=mean(dmL,dims=1)[1,:]
np=pyimport("numpy")

using LinearAlgebra

covDm=np["cov"](copy(transpose(dmL)))
eigVD,eVectD=eigen(covDm)
diag=Diagonal(1 ./eigVD);



invCov2=(eVectD)*diag*transpose(eVectD)



piaKu1,rrate1,zSim1=simZ(dmL[1,:])
radarMod.setZObs(copy(zSimL[1,:]))
radarMod.setInvCov(invCov2,dmMV)

function f_z(Dm)
    piaKu,rrate,zSim=simZ(Dm)
    n=size(zSim)[1]
    fobj=0.0
    for i=1:n
        fobj=fobj+(zSim[i]-Z_Obs[i])^2
    end
    #println(fobj)
    for i=1:n
        for j=1:n
            fobj=fobj+0.1*(Dm[i]-xMean[i])*invCov[i,j]*(Dm[j]-xMean[j])
        end
    end
    return fobj
end

function g_z(gradZ_out,Dm)
    piaKu,rrate,zSim,gradZ=simZ_g(Dm)
    n=size(zSim)[1]
    fobj=0.0
    for i=1:n
        fobj=fobj+(zSim[i]-Z_Obs[i])^2
    end
    gradZ_out.=2*gradZ*(zSim-Z_Obs)
    gradZ_out.=gradZ_out .+ 0.2*invCov*(Dm-xMean)
    for i=1:n
        for j=1:n
            fobj=fobj+0.1*(Dm[i]-xMean[i])*invCov[i,j]*(Dm[j]-xMean[j])
        end
    end
    return gradZ,zSim,piaKu,rrate
    #return fobj
end

fValue=f_z(dmL[1,:])
using Optim

using PyPlot
plot(Z_Obs,h1[100:165])
#
xsol=dmMV
gradZ_out=zeros(66,66);
for it=1:6
    global xsol
    println(f_z(xsol))
    gradZ,zSim,piaKu,rrate=g_z(gradZ_out,xsol)
    #gradZ=transpose(gradZ);
    dy=transpose(gradZ)*(Z_Obs-zSim)-invCov*(xsol-xMean);
    A=transpose(gradZ)*gradZ+invCov;
    dx=A\dy;
    xsol=xsol+1*dx;
end
gradZ,zSimf,piaKuf,rretf=g_z(gradZ_out,xsol)

#sol=Optim.minimizer(res)
#res=optimize(f_z,g_z,dmMV,method=BFGS(),iterations=10)
#sol=Optim.minimizer(res)
#piaKu,rrateL,zSim=simZ(sol);
#plot(zSim,h1[100:165])

