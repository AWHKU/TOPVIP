cd /home/athchu/data/arDCA
conda activate covmut

julia -t 8
mypkgdir = normpath(joinpath(pwd(),".."))
datadir=joinpath(mypkgdir,"data") # put here your path
using Pkg
Pkg.activate(mypkgdir)
using ArDCA


fasta_file="data/zeroshot_prediction/EVcoupling_alignments/UBE"
arnet,arvar=ardca(fasta_file, verbose=true, lambdaJ=0.02,lambdaH=0.001);

#θ = 0.3767497842212073 threshold = 349.0
#M = 52 N = 928 Meff = 45.0

target_sequence = 1

using ExtractMacro: @extract
@extract arnet:H J p0 idxperm
### H=76-element Vector{Vector{Float64}}
### J=76-element Vector{Array{Float64, 3}} 
### p0=21-element Vector
### idxperm 76-element Vector{Int64}
@extract arvar:Z M N q
### Z= 76×6152 Matrix, that is the SpCas9 length and the No. of seqs in the alignment
### M=6152
### N=76
### q=21


pc=0.1
ppc = (1 - pc) * p0 + pc * ones(q) / q
### 21-element Vector{Float64}

Da = fill(Inf64, q, N)
### 21×76 Matrix{Float64}

seqid=1
xori = Z[:, seqid]
### 928 element 
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
### 84.91360344290442 , WT arlike0

### So inorder to performa multimutation.
### first figure out which opsition=1318 and the AA 
### look at idxperm? - no. are not following aa pos because the column are sorted by conservation aka entropic order

### make alphebet ###
nums= [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
letters=[ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"]
LN=Dict()
LN=Dict( letters .=> nums )


using Pandas
df = read_csv("data/zeroshot_prediction/arDCA/UBE_arDCA_input.csv")
AA=df[:converted_mutations] # get AACombo
E=zeros(length(AA))


### for each mutation, split them, identify index, find mutant, find idxperm, swap 
function extract_muts(g)
	g=split(g, "_")
	mut_aa=Vector{Any}(nothing, length(g))
	pos=Vector{Any}(nothing, length(g))
	aa_index=Vector{Any}(nothing, length(g))
	for a in 1:length(g)
		i=g[a]
		k=string(i[length(i)])
		mut_aa[a]=k
		p=parse(Int64, i[2:length(i)-1])
		pos[a]=p
		idx=findall(x -> x == p, idxperm)
		aa_index[a]=idx
	end
	return mut_aa, pos, aa_index
end




for w in 1:length(AA)
	g=AA[w]
	if g !=NaN && g != "WT"
		mut_aa, pos, aa_index=extract_muts(g)
		arlike = zeros(N)
		xmut = copy(xori) ### rewrite xmut
		for i in 1:length(pos)
			AA_v=mut_aa[i]
			xmut[aa_index[i]] .= LN[AA_v]
		end
		_outputarnet!(arlike, xmut, J, H, ppc, N, q)
		E[w]=-sum(log.(arlike)) - ll0
	else
		E[w]=0
	end 
end

output = DataFrame(Dict(:AACombo=>AA, :arDCA=>E))
to_csv(output, "UBE_arDCA.csv")

