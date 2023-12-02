cd /home/athchu/data/arDCA
conda activate covmut

julia -t 8
mypkgdir = normpath(joinpath(pwd(),".."))
datadir=joinpath(mypkgdir,"data") # put here your path
using Pkg
Pkg.activate(mypkgdir)
using ArDCA


fasta_file="data/zeroshot_prediction/EVcoupling_alignments/SpCas9/SPCAS9_e40.a2m"
arnet,arvar=ardca(fasta_file, verbose=true, lambdaJ=0.02,lambdaH=0.001);

#function dms_single_site

target_sequence = 1
D,idxgap=dms_single_site(arnet,arvar,target_sequence)


using ExtractMacro: @extract
@extract arnet:H J p0 idxperm
### H=1317-element Vector{Vector{Float64}}
### J=1317-element Vector{Array{Float64, 3}} 
### p0=21-element Vector
### idxperm 1318-element Vector{Int64}
@extract arvar:Z M N q
### Z= 1318×679 Matrix, that is the SpCas9 length and the No. of seqs in the alignment
### M=679
### N=1318
### q=21

using DelimitedFiles
pc=0.1
ppc = (1 - pc) * p0 + pc * ones(q) / q
### 21-element Vector{Float64}

Da = fill(Inf64, q, N)
### 21×1318 Matrix{Float64}

seqid=1
xori = Z[:, seqid]
### 1318 element 
### that the encoding of the WT sequence
```
  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```

xmut = copy(xori)
###

idxnogap = findall(x -> x != 21, xori)
### remove location that has a gap in the WT sequence

arlike = zeros(N)
arlike0 = zeros(N)

using LoopVectorization: @avx 
DtotH = Dict{Tuple{Int,Int},Vector{Float64}}()

function softmax!(r::Vector{Float64}, x::Vector{Float64})
    n = length(x)
    u = maximum(x)
    s = 0.0
    @inbounds @avx for i = 1:n
        r[i] = exp(x[i] - u)
        s += r[i]
    end
    invs = inv(s)
    @inbounds @avx for i = 1:n
        r[i] *= invs
    end
    r
end

softmax!(x::Vector{Float64}) = softmax!(x, x)
softmax(x::Vector{Float64}) = softmax!(similar(x, Float64), x)



global _outputarnet!

function _outputarnet!(dest, x, J, H, p0, N, q)
        dest[1] = p0[x[1]]
        totH = Base.get!(DtotH,(q,Threads.threadid()),Vector{Float64}(undef,q))
        #fill!(totH,0.0)
        @inbounds for site in 1:N-1
            Js = J[site]
            h = H[site]
            copy!(totH,h)
            @avx for i in 1:site
                for a in 1:q
                    totH[a] += Js[a,x[i],i]
                end
            end
            softmax!(totH)
            dest[site+1]=totH[x[site+1]]
        end
        dest
end



_outputarnet!(arlike0, xori, J, H, ppc, N, q)

ll0 = -sum(log.(arlike0))
### 724.4040298398668 , WT arlike0

### make single mutants for each position

@inbounds for i in idxnogap
    if xori[i] == 21
        continue
    end
    for a in 1:q
    	if a != xori[i] # for a mutation
        	xmut[i] = a ### substitute mutation ###
        	_outputarnet!(arlike, xmut, J, H, ppc, N, q)
            Da[a, i] = -sum(log.(arlike)) - ll0 ## -log(mut) - log (WT)
        else
            Da[a, i] = 0.0
        end
    end
    xmut[i] = xori[i] #reset xmut to the original velue 
end
invperm = sortperm(idxperm)
D=Da[:, invperm]
idxgap=sort!(idxperm[setdiff(1:N, idxnogap)])


### So inorder to performa multimutation.
### first figure out which opsition=1318 and the AA 
### look at idxperm? - no. are not following aa pos because the column are sorted by conservation aka entropic order

### go to EVcoupling .align_standard.outcfg to find out the no to the corresponding to AApos
findall(x -> x == 644, idxperm) 
### 218
findall(x -> x == 675, idxperm)
### 263
findall(x -> x == 827, idxperm)
### 246
findall(x -> x == 902, idxperm)
### 172
findall(x -> x == 903, idxperm)
### 90
findall(x -> x == 905, idxperm)
### 71
findall(x -> x == 974, idxperm)
### 896


using Pandas
df = read_csv("data/zeroshot_prediction/EvoEF/SpCas9_R661Q695K848E923T924Q926K1003R1060_EvoEF.csv")
AA=df[:AACombo] # get AACombo


### make alphebet ###
nums= [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
letters=[ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"]
AA_index=[218, 263, 246, 172, 90, 71, 896]
LN=Dict()
LN=Dict( letters .=> nums )

E=zeros(length(AA))
for w in 1:length(AA)
	AA_v=split(AA[w], "")
	arlike = zeros(N)
	xmut = copy(xori) ### rewrite xmut
	for i in 1:7
		xmut[AA_index[i]] = LN[AA_v[i]]
	end
	_outputarnet!(arlike, xmut, J, H, ppc, N, q)
	E[w]=-sum(log.(arlike)) - ll0 
end

output = DataFrame(Dict(:AACombo=>AA, :arDCA=>E))
to_csv(output, "SpCas9_R661Q695K848E923T924Q926K1003R1060_arDCA.csv")


### SpG
findall(x -> x == 1135, idxperm) 
### 631
findall(x -> x == 1136, idxperm)
### 1191
findall(x -> x == 1218, idxperm)
### 780
findall(x -> x == 1219, idxperm)
### 320
findall(x -> x == 1335, idxperm)
### NA
findall(x -> x == 1337, idxperm)
### NA


using Pandas
df= read_csv("SpG_D1135S1136G1218E1219R1335T1337_EVmutation_ZeroShotPred.csv")
AA=df[:Combo] # get AACombo


### make alphebet ###
nums= [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
letters=[ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"]
AA_index=[631, 1191, 780, 320]
LN=Dict()
LN=Dict( letters .=> nums )

E=zeros(length(AA))
for w in 1:length(AA)
    AA_v=split(AA[w], "")
    arlike = zeros(N)
    xmut = copy(xori) ### rewrite xmut
    for i in 1:length(AA_index)
        xmut[AA_index[i]] = LN[AA_v[i]]
    end
    _outputarnet!(arlike, xmut, J, H, ppc, N, q)
    E[w]=-sum(log.(arlike)) - ll0 
end

output = DataFrame(Dict(:AACombo=>AA, :arDCA=>E))
to_csv(output, "SpG_D1135S1136G1218E1219R1335T1337_arDCA.csv")

