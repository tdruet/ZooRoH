! PROGRAM : ZooRoH
! Author  : Tom DRUET
! Copyright (C) 2017

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!See <http://www.gnu.org/licenses/>.

program ZooRoH
implicit none
integer ::i,j,k,l,io,n,npos,nclust,round,nind,id,maxid,bp,chr,testid,testid2,nF,nl
integer ::lf,nhet,nsnps,starts,ends,nchr,lpos,fpos,niter,ift
integer ::readf,estimateG,onerate,nparam
integer*2, allocatable ::geno(:),genos(:,:)
integer, allocatable ::posi(:),chr_limits(:,:),isF(:),IT(:,:)
real*8, allocatable ::alpha(:,:),scaling(:),beta(:,:),gamma(:,:),ems(:,:),freq(:),pinit(:),Fs(:),hsprobs(:),as(:),ainit(:),Finit(:)
real*8, allocatable ::outprob(:,:,:),num_bjk(:,:),den_bjk(:),num_pi(:),den_pi(:),trans(:,:),T2(:,:),homoz(:),freqin(:),pemission(:,:)
real*8, allocatable ::globalF(:,:)
real*8 ::f0,f1,f2,ptrans,lik,val,num_transition,pos,maxprob,pF,perr,loglik,maf,minmaf,AIC,BIC
real*8 ::tmp2,lseg,fromitoj,minrate,maxrate
real*8, allocatable ::global(:),fromi(:),ini(:),mux(:),EL(:)
real*4, allocatable ::genosr(:,:),genor(:)
real*8, parameter ::Morgan=100000000.d0
character*50 ::mname,chrom,all1,all2,chromold,filename,freqfile

character*50 ::simfile,outfile,str1,str4,OUTPUT
integer ::checksim,nclass,maxclass
integer, allocatable ::classes(:,:),oldnum(:),statein(:)
real*8, allocatable :: checktable(:,:),probtable(:,:,:)

! file with genotypes at each marker for one individual: number of 1 alleles (0 for 22, 1 for 12 and 2 for 11)
! frequencies: frequency of allele 1

!############ INITIALIZE PARAMETERS ################

call read_parameters

open(10,file=filename) 
open(11,file='segmentsF.txt')
open(12,file='mixingF.txt')
open(13,file='averageF.txt')
open(14,file='countsF.txt')
if(estimateG==1)open(15,file='ageF.txt')
!open(16,file='inbreeding.prob')

npos=0;chromold='chr0';nchr=0
do
read(10,*,iostat=io)chrom,mname,bp
if(io/=0)exit
npos=npos+1
if(chrom .ne. chromold)nchr=nchr+1
chromold=chrom
enddo

if(ift==1)allocate(geno(nind),freq(npos),posi(npos),genos(nind,npos),chr_limits(nchr,2))
if(ift==2 .or. ift==3)allocate(genor(3*nind),freq(npos),posi(npos),genosr(3*nind,npos),chr_limits(nchr,2))
if(ift==4)allocate(geno(2*nind),freq(npos),posi(npos),genos(2*nind,npos),chr_limits(nchr,2))
allocate(homoz(nind))
if(OUTPUT .eq. "PartialF")allocate(outprob(npos,nind,sum(isF)))
if(OUTPUT .eq. "PartialF" .or. OUTPUT .eq. "TotalF")allocate(globalF(npos,nind))
if(OUTPUT .eq. "PartialF")outprob=0.00
if(OUTPUT .ne. "no")globalF=0.00

if(checksim==1)then
 allocate(classes(nind,npos),oldnum(npos),statein(nind))
 maxclass=0;classes=0;oldnum=0
 open(101,file=simfile)
 k=0
 do 
  read(101,*,iostat=io)chrom,bp,all1,all2,(statein(i),i=1,nind)
  if(io/=0)exit
  k=k+1
  classes(:,k)=statein
  maxclass=max(maxclass,maxval(statein))
 enddo
 close(101)
 print*,'Number of simulated classes ::',maxclass
 allocate(checktable(maxclass,nclust),probtable(100,nclust,2))
 checktable=0.d0;probtable=0.d0
 deallocate(statein)
endif

rewind(10)

if(readf==1)then
 allocate(freqin(npos))
 freqin=0.d0
 k=0
 open(101,file=freqfile)
 do
  read(101,*,iostat=io)f0
  if(io/=0)exit
  k=k+1
  if(k > npos)then
    print*,'Number of marker in frequence files larger than in genotype file !!!'
    exit
  endif
  freqin(k)=f0
 enddo
 close(101)
endif

n=0;k=0;posi=0.00;chr_limits=0;chromold='chr0';chr=0;freq=0.d0
do
 if(ift==1)read(10,*,iostat=io)chrom,mname,bp,all1,all2,(geno(i),i=1,nind)
 if(ift==2 .or. ift==3)read(10,*,iostat=io)chrom,mname,bp,all1,all2,(genor(i),i=1,3*nind)
 if(ift==4)read(10,*,iostat=io)chrom,mname,bp,all1,all2,(geno(i),i=1,2*nind)
 if(io/=0)exit
 k=k+1
 if(readf==0)then
  call get_freq
 else
  f0=freqin(k)
 endif
 maf=min(f0,1.d0-f0)
 if(maf<minmaf)cycle
 n=n+1
 if(checksim==1)oldnum(n)=k
! posi(n)=1.d0*bp/1000000.0
 posi(n)=bp
 if(ift==1 .or. ift==4)genos(:,n)=geno
 if(ift==2 .or. ift==3)genosr(:,n)=genor
 freq(n)=f0
 if(chrom .ne. chromold)then
   chr=chr+1
   chr_limits(chr,1)=n
   chromold=chrom
 endif
 chr_limits(chr,2)=n
enddo
print*,'Number of variants ::',npos
print*,'Number of selected variants (based on MAF) ::',n
npos=n

if(readf==1)deallocate(freqin)

call esthomoz

if(ift==1)deallocate(geno)
if(ift==1)allocate(geno(npos))

print*,'Number of individuals red ::',nind

allocate(trans(nclust,nclust),pinit(nclust),T2(nclust,nclust),num_pi(nclust),Fs(nclust),hsprobs(nF),pemission(npos,2))
if(estimateG==1)allocate(global(nclust),fromi(nclust),ini(nclust),EL(nclust),mux(nclust))

!###### start HMM

allocate(alpha(nclust,npos),scaling(npos))
allocate(beta(nclust,npos),gamma(nclust,npos))

do id=testid,testid2
 Fs=Finit;as=ainit
 if(ift==1)geno=genos(id,:)
 call get_pemission(id)
 do round=1,niter
 pinit=Fs


alpha=0.d0;scaling=0.d0
beta=0.d0;gamma=0.d0
num_pi=0.d0
if(estimateG==1)then
 global=0.d0;EL=0.d0
endif
loglik=0.d0

!############ FORWARD ALGORITHM ####################

do chr=1,nchr
fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! initialisation

do i=1,nclust
  alpha(i,fpos)=pinit(i)*pemission(fpos,isF(i)+1)
  scaling(fpos)=scaling(fpos)+alpha(i,fpos)
enddo

 scaling(fpos)=1.0/scaling(fpos)
 alpha(:,fpos)=alpha(:,fpos)*scaling(fpos)

! induction: to get to i, two ways:
! 1) transition + jump into cluster i
! 2) no transition and in i at previous position
 
 do k=fpos+1,lpos
  trans=transition(k-1)
  do i=1,nclust
   do j=1,nclust
    alpha(i,k)=alpha(i,k)+alpha(j,k-1)*trans(j,i)*pemission(k,isF(i)+1)
   enddo
    scaling(k)=scaling(k)+alpha(i,k)
  enddo
 scaling(k)=1.0/scaling(k)
 alpha(:,k)=alpha(:,k)*scaling(k)
 enddo  

! termination

 loglik=loglik-sum(log(scaling(fpos:lpos)))

!############ BACKWARD ALGORITHM ####################


! initialisation

do i=1,nclust
 gamma(i,lpos)=alpha(i,lpos)*1.0  ! beta(i,j,npos)=1.0
 beta(i,lpos)=1.0*scaling(lpos)
enddo

! induction
! to arrive in k: with or without transition

do k=lpos-1,fpos,-1
 trans=transition(k)
 do i=1,nclust
  do j=1,nclust
    beta(i,k)=beta(i,k)+trans(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
  enddo
  beta(i,k)=beta(i,k)*scaling(k)
  gamma(i,k)=alpha(i,k)*beta(i,k)/scaling(k)
 enddo
enddo

do k=fpos,lpos
 gamma(:,k)=gamma(:,k)/sum(gamma(:,k))
enddo

!########## estimation transitions #############

do i=1,nclust
 num_pi(i)=num_pi(i)+gamma(i,fpos)
enddo

if(estimateG==0)then
 do k=fpos,lpos-1
  T2=pexit(k)
  do i=1,nclust
   do j=1,nclust
    num_pi(j)=num_pi(j)+alpha(i,k)*T2(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
   enddo
  enddo
 enddo
else
 do k=fpos,lpos-1
  T2=pexit(k)
  trans=transition(k)
  fromi=0.d0;ini=0.d0;tmp2=0.d0
  do i=1,nclust
   lseg=(posi(k+1)-posi(k))/Morgan
   if(lseg>0.d0)mux(i)=1/as(i)-lseg/(dexp(as(i)*lseg)-1.d0)
   do j=1,nclust
    fromitoj=alpha(i,k)*T2(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
    num_pi(j)=num_pi(j)+fromitoj
    ini(i)=ini(i)+alpha(i,k)*trans(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
!    tmp2=tmp2+alpha(i,k)*trans(i,j)*emission(j,k+1)*beta(j,k+1)
    fromi(i)=fromi(i)+fromitoj
   enddo
  enddo
  do i=1,nclust
   global(i)=global(i)+fromi(i)
!   if(lseg>0.d0)EL(i)=EL(i)+((ini(i)-fromi(i))/tmp2)*(1.d0*(posi(k+1)-posi(k))/100000000.d0)+(fromi(i)/tmp2)*mux(i)
   if(lseg>0.d0)EL(i)=EL(i)+(ini(i)-fromi(i))*(1.d0*(posi(k+1)-posi(k))/Morgan)+fromi(i)*mux(i)
  enddo
 enddo
endif

!print'(i6,1x,i2,<nclust>(1x,f19.7))',id,round,num_pi

enddo ! end loop on chromosomes

if(estimateG==1)then
 if(nclust==2 .and. onerate==1)then ! joint parameter
   as(1)=sum(global)/sum(EL)
   if(as(1)<minrate)as(1)=minrate
   if(as(1)>maxrate)as(1)=maxrate
   as(2)=as(1)
 else
  do i=1,nclust
    as(i)=global(i)/EL(i)
   if(as(i)<minrate)as(i)=minrate
   if(as(i)>maxrate)as(i)=maxrate
  enddo
 endif
endif

!if(mod(round,1)==0)print'(i5,1x,f15.6,*(1x,f15.6))',round,loglik,(as(i),i=1,nclust)


Fs=num_pi/sum(num_pi)
if(round==niter)then
 if(estimateG==1)nparam=2*nclust-1
 if(estimateG/=1)nparam=nclust-1
 if(estimateG==1 .and. onerate==1 .and. nclust==2)nparam=2
 AIC=-2.d0*loglik+2.d0*nparam
 BIC=-2.d0*loglik+log(1.d0*npos)*nparam
 print'(i6,1x,i4,3(1x,f15.6),*(1x,f9.7,1x,f11.4))',id,round,loglik,AIC,BIC,(Fs(i),as(i),i=1,nclust)
 if(estimateG==1)write(15,'(i5,*(1x,f11.4))')id,(as(i),i=1,nclust)
endif

enddo

l=0
do k=1,nclust
 if(isF(k)==0)cycle
 l=l+1
 if(OUTPUT .eq. 'PartialF')outprob(1:npos,id,l)=gamma(k,:)
 if(OUTPUT .ne. 'no')globalF(1:npos,id)=globalF(1:npos,id)+gamma(k,:)
 hsprobs(l)=sum(gamma(k,:))/(1.d0*npos)
enddo

write(12,'(i6,*(1x,f9.6))')id,(Fs(j),j=1,nclust)
write(13,'(i6,*(1x,f9.6))')id,(hsprobs(j),j=1,nF),sum(hsprobs),homoz(id)
write(14,'(i6,*(1x,f9.2))')id,(num_pi(j),j=1,nclust)


do chr=1,nchr
starts=0;ends=0;maxprob=0.00
do k=chr_limits(chr,1),chr_limits(chr,2)
! write(16,'(i6,1x,i2,1x,i9,1x,f6.4,*(1x,f9.6))')id,chr,posi(k),freq(k),(gamma(l,k),l=1,nclust)
 pF=dot_product(gamma(:,k),isF)
 if(pF > 0.99)then
  if(starts==0)starts=k
  maxprob=max(maxprob,pF)
  ends=k
 else if(pF <= 0.99)then
  if(starts /= 0 .and. maxprob > 0.999)then
   nhet=0;nsnps=0
   do l=starts,ends
    if(ift==1 .and. geno(l)==1)nhet=nhet+1
    if(ift==1 .and. geno(l)/=9)nsnps=nsnps+1
   enddo
!   nsnps=ends-starts+1
   lf=posi(ends)-posi(starts)+1
   nl=0
   do l=1,nclust
     if(isF(l)==0)cycle
     nl=nl+1
     hsprobs(nl)=sum(gamma(l,starts:ends))
   enddo
   hsprobs=hsprobs/sum(hsprobs)
!   write(11,'(i4,1x,i2,2(1x,i7),3(1x,i9),2(1x,i7),<nF+1>(1x,f9.6))')id,chr,starts,ends,posi(starts),posi(ends),lf,nhet,nsnps,maxprob,(hsprobs(l),l=1,nF)
   write(11,'(i4,1x,i2,2(1x,i7),3(1x,i9),2(1x,i7),*(1x,f9.6))')id,chr,starts,ends,posi(starts),posi(ends),lf,nhet,nsnps,maxprob,(hsprobs(l),l=1,nF)
  endif
  starts=0;ends=0;maxprob=0.00
 endif  
enddo

if(starts /=0 .and. maxprob > 0.999)then
 nhet=0;nsnps=0
 do l=starts,ends
  if(ift==1 .and. geno(l)==1)nhet=nhet+1
  if(ift==1 .and. geno(l)/=9)nsnps=nsnps+1
 enddo
! nsnps=ends-starts+1
 lf=posi(ends)-posi(starts)+1
 nl=0
 do l=1,nclust
   if(isF(l)==0)cycle
   nl=nl+1
   hsprobs(nl)=sum(gamma(l,starts:ends))
 enddo
 hsprobs=hsprobs/sum(hsprobs)
! write(11,'(i4,1x,i2,2(1x,i7),3(1x,i9),2(1x,i7),<nF+1>(1x,f9.6))')id,chr,starts,ends,posi(starts),posi(ends),lf,nhet,nsnps,maxprob,(hsprobs(l),l=1,nF)
 write(11,'(i4,1x,i2,2(1x,i7),3(1x,i9),2(1x,i7),*(1x,f9.6))')id,chr,starts,ends,posi(starts),posi(ends),lf,nhet,nsnps,maxprob,(hsprobs(l),l=1,nF)
endif
enddo ! output per chr

if(checksim==1)then
 do k=1,npos
  do i=1,nclust
    checktable(classes(id,oldnum(k)),i)=checktable(classes(id,oldnum(k)),i)+gamma(i,k)
    j=ceiling(gamma(i,k)*100.d0)
    if(j<1)j=1
    if(j>100)j=100
    if(classes(id,oldnum(k))==1)probtable(j,i,1)=probtable(j,i,1)+1.d0
    probtable(j,i,2)=probtable(j,i,2)+1.d0
  enddo
 enddo
endif

enddo ! id

if(checksim==1)then
 open(101,file='sim_by_res.txt')
 open(102,file='sim_by_probres.txt')
 do i=1,nclust
   write(101,'(i2,*(1x,f14.4))')i,(checktable(j,i),j=1,maxclass)
 enddo
 do j=1,100
   do i=1,nclust
     if(probtable(j,i,2)>0.d0)probtable(j,i,1)=probtable(j,i,1)/probtable(j,i,2)
   enddo
   write(102,'(i3,*(1x,f7.5,1x,f14.4))')j,(probtable(j,i,1),probtable(j,i,2),i=1,nclust)
 enddo
endif


! ****** OUTPUT *********

if(OUTPUT .eq. 'PartialF')then
l=0
do i=1,nclust
 if(isF(i)==1)then
  l=l+1
  if(l<10)write(str4,'(i1)')l
  if(l>9 .and. l<100)write(str4,'(i2)')l
  if(l>99 .and. l<1000)write(str4,'(i3)')l
  str1="PartialF"//trim(str4)//".txt"
  outfile=str1//".txt"
  open(112,file=outfile)
  do chr=1,nchr
   do k=chr_limits(chr,1),chr_limits(chr,2)
    write(112,'(i8,1x,i10,1x,i3,*(1x,f8.6))')k,posi(k),chr,(outprob(k,id,l),id=testid,testid2)
   enddo
  enddo
  close(112)
 endif
enddo
endif

if(OUTPUT .ne. 'no')then
  open(112,file='TotalF.txt')
  do chr=1,nchr
   do k=chr_limits(chr,1),chr_limits(chr,2)
    write(112,'(i8,1x,i10,1x,i3,*(1x,f8.6))')k,posi(k),chr,(globalF(k,id),id=testid,testid2)
   enddo
  enddo
  close(112)
endif

contains


function transition(marker)
implicit none
integer ::marker,hs,hs2
real*8 ::transition(nclust,nclust),sumF(2),F,a,r,gr

!as=(/2.00,2.00,1.00,0.50,0.10,0.05/)
!as=(/2.00,200.00,20.00,2.00,0.20,0.05/)

!sumF(1)=sum(Fs,isF == 0)
!sumF(2)=sum(Fs, isF == 1)

transition=0.d0
do hs=1,nclust
!gr=isF(hs)+1
sumF(1)=dot_product(Fs,IT(hs,:))
a=as(hs)
r=dexp(-a*(posi(marker+1)-posi(marker))/Morgan)
 do hs2=1,nclust
!   if(isF(hs)==isF(hs2))cycle ! only from background to F and vice versa
!   transition(hs,hs2)=(1-r)*Fs(hs2)/sumF(gr)
   if(IT(hs,hs2)==0)cycle
   transition(hs,hs2)=(1-r)*Fs(Hs2)/sumF(1)
 enddo
 transition(hs,hs)=transition(hs,hs)+r
enddo

end function

function pexit(marker)
implicit none
integer ::marker,hs,hs2
real*8 ::pexit(nclust,nclust),sumF(2),a,r,gr

!as=(/2.00,2.00,1.00,0.50,0.10,0.05/)
!as=(/2.00,200.00,20.00,2.00,0.20,0.05/)

!sumF(1)=sum(Fs,isF == 0)
!sumF(2)=sum(Fs, isF == 1)

pexit=0.d0
do hs=1,nclust
!gr=isF(hs)+1
sumF(1)=dot_product(Fs,IT(hs,:))
a=as(hs)
r=dexp(-a*(posi(marker+1)-posi(marker))/Morgan)
 do hs2=1,nclust
!   if(isF(hs)==isF(hs2))cycle
!   pexit(hs,hs2)=(1-r)*Fs(hs2)/sumF(gr)
   if(IT(hs,hs2)==0)cycle
   pexit(hs,hs2)=(1-r)*Fs(Hs2)/sumF(1)
 enddo
enddo

end function

subroutine read_parameters
implicit none
integer :: check
character*200 ::inputline,label1,iformat

open(20,file='param.txt')

nclust=0;testid=0;testid2=0;nind=0;ift=1
minmaf=0.d0;readf=0;estimateG=0;niter=100;onerate=0
OUTPUT='no'

print*,'        '
print*,'/***** PARAMETER FILE *****/'
print*,'        '

do
read(20,*,iostat=io)inputline
if(io/=0)exit
check=0
print*,inputline
if(inputline(1:18) .eq. "#NUMBER_OF_CLASSES")then
 read(20,*)nclust
 print*,nclust
 allocate(as(nclust),ainit(nclust),isF(nclust),Finit(nclust),IT(nclust,nclust))
 IT=1;as=0.d0;isF=0;Finit=0.d0
 check=1
endif
if(inputline(1:16) .eq. "#RATE_PARAMETERS")then
 if(nclust==0)then
   print*,'NUMBER OF CLASSES shoud come first !'
   stop
 endif
 read(20,*,iostat=io)(ainit(i),i=1,nclust)
 print'(*(1x,f9.6))',(ainit(i),i=1,nclust)
 check=1
endif
if(inputline(1:22) .eq. "#INBREEDING_INDICATORS")then
 if(nclust==0)then
   print*,'NUMBER OF CLASSES shoud come first !'
   stop
 endif
 read(20,*,iostat=io)(isF(i),i=1,nclust)
 nF=sum(isF)
 print'(*(1x,i1))',(isF(i),i=1,nclust)
 check=1
endif
if(inputline(1:28) .eq. "#STARTING_MIXING_PROPORTIONS")then
 if(nclust==0)then
   print*,'NUMBER OF CLASSES shoud come first !'
   stop
 endif
 read(20,*,iostat=io)(Finit(i),i=1,nclust)
 print'(*(1x,f9.6))',(Finit(i),i=1,nclust)
 check=1
endif
if(inputline(1:18) .eq. "#TRANSITION_MATRIX")then
 if(nclust==0)then
   print*,'NUMBER OF CLASSES shoud come first !'
   stop
 endif
 do k=1,nclust
  read(20,*,iostat=io)(IT(k,i),i=1,nclust)
  if(io/=0)exit
  print'(*(1x,i1))',(IT(k,i),i=1,nclust)
 enddo
 check=1
endif
if(inputline(1:11) .eq. "#ERROR_RATE")then
 read(20,*,iostat=io)perr
 print*,perr
 check=1
endif
if(inputline(1:22) .eq. "#NUMBER_OF_INDIVIDUALS")then
 read(20,*,iostat=io)nind
 print*,nind
 check=1
endif
if(inputline(1:11) .eq. "#INPUT_FILE")then
 read(20,*,iostat=io)filename
 print*,filename
 check=1
endif
if(inputline(1:15) .eq. "#ANALYSIS_RANGE")then
 read(20,*,iostat=io)testid,testid2
 print*,testid,testid2
 check=1
endif
if(inputline .eq. "#MINMAF")then
  read(20,*,iostat=io)minmaf
  print*,minmaf
  check=1
endif
if(inputline .eq. "#FREQUENCIES")then
  read(20,*,iostat=io)freqfile
  print*,freqfile
  readf=1
  check=1
endif
if(inputline .eq. "#ITERATIONS")then
  read(20,*,iostat=io)niter
  print*,niter
  check=1
endif
if(inputline .eq. "#ESTIMATE_RATES")then
  read(20,*,iostat=io)minrate,maxrate
  estimateG=1
  check=1
endif
if(inputline .eq. "#ONE_RATE")then
  onerate=1
  check=1
endif
if(inputline .eq. "#SIMULATED_CLASSES")then
  checksim=1
  read(20,*,iostat=io)simfile
  print*,simfile
  check=1
endif
if(inputline .eq. "#INPUT_FORMAT")then
  read(20,*,iostat=io)iformat
  if(iformat=='GEN')ift=1
  if(iformat=='GP')ift=2
  if(iformat=='GL')ift=3 !## three floating point?
  if(iformat=='AD')ift=4 !## allele depth (per allele)
  print*,iformat,ift
  if((iformat .ne. 'GEN') .and. (iformat .ne. 'GP') .and. (iformat .ne. 'GL') .and. (iformat .ne. 'AD'))then
    print*,'Unknown input format ::',iformat
    print*,'Should ne GEN, GP, GL or AD'
    print*,'Default format will be used'
  endif
  check=1  
endif
if(inputline .eq. "#OUTPUT")then
  read(20,*,iostat=io)OUTPUT
  print*,OUTPUT
  if(OUTPUT .ne. 'PartialF' .and. OUTPUT .ne. 'TotalF')then
    print*,'Warning: output format should be PartialF or TotalF !'
    print*,'Format ignored !'
    OUTPUT='no'
  endif
  check=1
endif
if(check==0)print*,'Unknown option (ignored) ::',inputline
if(io/=0)exit
enddo

if(testid<=0)testid=1
if(testid2<=0 .or. testid2>nind)testid2=nind

print*,''
print*,'/***** end parameter file *****/'
print*,''

end subroutine

subroutine get_pemission(id)
integer ::id,k
real*8 ::maf,f1,f2,pr11,pr12,pr22,prT

do k=1,npos
 pemission(k,:)=1.d0
 
 f1=freq(k);f2=1.0-freq(k)
 maf=min(f1,f2)

 pr11=0.d0;pr12=0.d0;pr22=0.d00
 if(ift==1)then !### genos = number of alleles 1
  if(genos(id,k)==9)cycle
  if(genos(id,k)==0)pr22=1.d0
  if(genos(id,k)==1)pr12=1.d0
  if(genos(id,k)==2)pr11=1.d0
 endif

 if(ift==2)then !### genoprob = prob(11), prob(12) & prob(22)
  if((genosr(3*id-2,k)+genosr(3*id-1,k)+genosr(3*id,k))==0.d00)cycle
  pr11=genosr(3*id-2,k)
  pr12=genosr(3*id-1,k)
  pr22=genosr(3*id,k)
 endif

 if(ift==3)then !### GL = GL11, GL12 & GL22
  if((genosr(3*id-2,k)+genosr(3*id-1,k)+genosr(3*id,k))==0.d00)cycle
  pr11=10**(-genosr(3*id-2,k)/10.d0)
  pr12=10**(-genosr(3*id-1,k)/10.d0)
  pr22=10**(-genosr(3*id,k)/10.d0)
  prT=pr11+pr12+pr22
  pr11=pr11/prT;pr12=pr12/prT;pr22=pr22/prT
 endif

 if(ift==4)then !#### AD => D allele 1 / D allele 2
  if((genos(2*id-1,k)+genos(2*id,k))==0)cycle
  pr11=(1.d0**genos(2*id-1,k))*(0.d0**genos(2*id,k))
  pr12=0.5**(genos(2*id-1,k)+genos(2*id,k))
  pr22=(0.d0**genos(2*id-1,k))*(1.d0**genos(2*id,k))
  prT=pr11+pr12+pr22
  pr11=pr11/prT;pr12=pr12/prT;pr22=pr22/prT
 endif

 pemission(k,1)=pr11*f1**2+pr12*2*f1*f2+pr22*f2**2
 pemission(k,2)=(pr11*f1+pr22*f2)*(1.d0-perr)+pr12*perr
enddo

end subroutine

subroutine get_freq
real*8 ::pr11,pr12,pr22,prT,r0,r1

l=0;r0=0.d0;f0=0.d0

do i=1,nind
   
 pr11=0.d0;pr12=0.d0;pr22=0.d00
 if(ift==1)then !### genos = number of alleles 1
  if(geno(i)==9)cycle
  if(geno(i)==0)pr22=1.d0
  if(geno(i)==1)pr12=1.d0
  if(geno(i)==2)pr11=1.d0
 endif

 if(ift==2)then !### genoprob = prob(11), prob(12) & prob(22)
  if((genor(3*i-2)+genor(3*i-1)+genor(3*i))==0.d00)cycle
  pr11=genor(3*i-2)
  pr12=genor(3*i-1)
  pr22=genor(3*i)
 endif

 if(ift==3)then !### GL = GL11, GL12 & GL22
  if((genor(3*i-2)+genor(3*i-1)+genor(3*i))==0.d00)cycle
  pr11=10**(-genor(3*i-2)/10.d0)
  pr12=10**(-genor(3*i-1)/10.d0)
  pr22=10**(-genor(3*i)/10.d0)
  prT=pr11+pr12+pr22
  pr11=pr11/prT;pr12=pr12/prT;pr22=pr22/prT
 endif

 if(ift==4)then !#### AD => D allele 1 / D allele 2
  if((geno(2*i-1)+geno(2*i)==0))cycle
  pr11=(1.d0**geno(2*i-1))*(0.d0**geno(2*i))
  pr12=0.5**(geno(2*i-1)+geno(2*i))
  pr22=(0.d0**geno(2*i-1))*(1.d0**geno(2*i))
  prT=pr11+pr12+pr22
  pr11=pr11/prT;pr12=pr12/prT;pr22=pr22/prT
  r1=1.d0*(geno(2*i-1))/(1.d0*(geno(2*i-1)+geno(2*i)))
 endif

   
 f0=f0+2.00*pr11+1.00*pr12
 if(ift==4)r0=r0+r1
 l=l+1

enddo

if(ift<4)f0=f0/(2.d0*l)
if(ift==4)f0=r0/(1.d0*l)

end subroutine

subroutine esthomoz
real*8 ::pr11,pr12,pr22,prT

homoz=0.d0

do id=1,nind
 n=0
 do k=1,npos
 
  f1=freq(k);f2=1.0-freq(k)
  maf=min(f1,f2)

  pr11=0.d0;pr12=0.d0;pr22=0.d00
  if(ift==1)then !### genos = number of alleles 1
   if(genos(id,k)==9)cycle
   if(genos(id,k)==0)pr22=1.d0
   if(genos(id,k)==1)pr12=1.d0
   if(genos(id,k)==2)pr11=1.d0
  endif

  if(ift==2)then !### genoprob = prob(11), prob(12) & prob(22)
   if((genosr(3*id-2,k)+genosr(3*id-1,k)+genosr(3*id,k))==0.d00)cycle
   pr11=genosr(3*id-2,k)
   pr12=genosr(3*id-1,k)
   pr22=genosr(3*id,k)
  endif

  if(ift==3)then !### GL = GL11, GL12 & GL22
   if((genosr(3*id-2,k)+genosr(3*id-1,k)+genosr(3*id,k))==0.d00)cycle
   pr11=10**(-genosr(3*id-2,k)/10.d0)
   pr12=10**(-genosr(3*id-1,k)/10.d0)
   pr22=10**(-genosr(3*id,k)/10.d0)
   prT=pr11+pr12+pr22
   pr11=pr11/prT;pr12=pr12/prT;pr22=pr22/prT
  endif

  if(ift==4)then !#### AD => D allele 1 / D allele 2
   if((genos(2*id-1,k)+genos(2*id,k))==0)cycle
   pr11=(1.d0**genos(2*id-1,k))*(0.d0**genos(2*id,k))
   pr12=0.5**(genos(2*id-1,k)+genos(2*id,k))
   pr22=(0.d0**genos(2*id-1,k))*(1.d0**genos(2*id,k))
   prT=pr11+pr12+pr22
   pr11=pr11/prT;pr12=pr12/prT;pr22=pr22/prT
  endif

  n=n+1
  homoz(id)=homoz(id)+pr11+pr22
  
 enddo
 if(n>0)homoz(id)=homoz(id)/(1.d0*n)
enddo

end subroutine


end program
