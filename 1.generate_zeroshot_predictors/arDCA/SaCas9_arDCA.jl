cd /home/athchu/data/arDCA
conda activate covmut




julia -t 8
mypkgdir = normpath(joinpath(pwd(),".."))
datadir=joinpath(mypkgdir,"data") # put here your path
using Pkg
Pkg.activate(mypkgdir)
using ArDCA


fasta_file="data/zeroshot_prediction/EVcoupling_alignments/SaCas9/J7RUA5_b0.1.a2m"
arnet,arvar=ardca(fasta_file, verbose=true, lambdaJ=0.02,lambdaH=0.001);

function dms_single_site

target_sequence = 1
D,idxgap=dms_single_site(arnet,arvar,target_sequence)



using ExtractMacro: @extract
@extract arnet:H J p0 idxperm
### H=1033-element Vector{Vector{Float64}}
### J=1032-element Vector{Array{Float64, 3}} 
### p0=21-element Vector
### idxperm 1033-element Vector{Int64}
@extract arvar:Z M N q
### Z= 1033×2472  Matrix, that is the SpCas9 length and the No. of seqs in the alignment
### M=2472
### N=1033
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
### 787.8454230582206 , WT arlike0

### So inorder to performa multimutation.
### first figure out which opsition=1318 and the AA 
### look at idxperm? - no. are not following aa pos because the column are sorted by conservation aka entropic order

### go to EVcoupling /media/achu/新增磁碟區/MLDE/2021_11_09_new_MLDE/SaCas9_MLDE/SaCas9_EVcoupling_results/SaCas9/align/SaCas9_align.outcfg to find out the no to the corresponding to AApos
# SaCas9_L887_N888_A889_N985_N986_L988_L989_R991
#```
#  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#```

findall(x -> x == 878, idxperm) 
### 224 
findall(x -> x == 879, idxperm)
### 658
findall(x -> x == 880, idxperm)
### 766
findall(x -> x == 976, idxperm)
### 985
findall(x -> x == 977, idxperm)
### 1001
findall(x -> x == 979, idxperm)
### 729
findall(x -> x == 980, idxperm)
### 915
findall(x -> x == 982, idxperm)
### 933


using Pandas
df = read_csv("data/zeroshot_prediction/arDCA/Sa1296_Fitness.csv")
AA=df[:AACombo] # get AACombo


### make alphebet ###
nums= [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
letters=[ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"]
AA_index=[224, 658, 766, 985, 1001, 729, 915, 933]
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
to_csv(output, "SaCas9_WED-PI_arDCA.csv")

### N888_A889_L909
findall(x -> x == 879, idxperm)
### 658
findall(x -> x == 880, idxperm)
### 766
findall(x -> x == 900, idxperm)
### 456

df = read_csv("data/zeroshot_prediction/arDCA/Sa8000_Fitness.csv")
AA=df[:AACombo] # get AACombo


### make alphebet ###
nums= [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
letters=[ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"]
AA_index=[658, 766, 456]
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
to_csv(output, "SaCas9_N888_A889_L909_arDCA.csv")
