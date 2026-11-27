### A Pluto.jl notebook ###
# v1.1.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f8f60370-e96f-11f1-9a04-21ec278d93a8
begin 
	
	import Pkg 
	Pkg.activate(".")
end

# ╔═╡ a31e3058-eac7-4d6e-9e3e-8376aa9543bf
begin 
	using ScaledPolynomials 
	using LinearAlgebra
	using DelimitedFiles
	using Plots
	using PlutoUI
	using Interpolations
	using Statistics
	using Observables
end

# ╔═╡ 6ba449bc-ff95-4bb2-b023-bd072e56c319
using NativeFileDialog

# ╔═╡ 06272905-c414-4a00-9373-1e216fc08906
begin 
	note_dir = @__DIR__
	ethalons_dir = joinpath(note_dir , "ethalons")
	data_dir = joinpath(note_dir , "data")
	data_dir_ir = joinpath(data_dir , "IR")
	data_dir_nir = joinpath(data_dir , "NIR")
	data_dir_vis = joinpath(data_dir , "VIS")
	results_dir = joinpath(note_dir , "results")
	isdir(results_dir) || mkdir(results_dir)
	cd(note_dir)
end

# ╔═╡ 9681c7ff-7512-4d80-8eee-6af7e7597ad4
html"""
<style>
.plutoui-rangeslider {
		width : 100%;
}
</style>
"""

# ╔═╡ af55b1c3-8825-46b1-986b-e2af7e57907a
plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 4 , xlabel = "Wavelength, μm" , ylabel = "R");

# ╔═╡ e5152a42-16c6-4ee2-a3a9-479b417caba3
vis_nir_merged = Ref{Matrix{Float64}}()

# ╔═╡ ee742d64-e0d6-40ed-9733-cec0d59031dd
fully_merged = Ref{Matrix{Float64}}()

# ╔═╡ 66373b23-b865-4b9f-9042-722bf6b9288f
path = save_file(@__DIR__)

# ╔═╡ ed3276ba-7c3e-42ea-88c6-4d58fe9c1d3a
p_merged = Ref(plot(; plot_common_args...));

# ╔═╡ aec51390-534c-4493-b1f8-e139032acbc9
begin 
	refresh_full_plot = Observable(false)
	on(refresh_full_plot) do _
		isassigned(fully_merged) || !isempty(fully_merged[]) || return nothing
		p_merged[] = plot(fully_merged[][:,1] ,  fully_merged[][:,2] ; plot_common_args...)
		p_merged[]
	end
end

# ╔═╡ 3cf28b79-09e3-4bed-9b9c-2811fbcdf4d7

function calibrate(V , eth)
	_V = copy(V)
	y = @view _V[:,2]
	x = @view _V[:,1]
	@. y *= eth(x)
	return _V
end

# ╔═╡ ab6c19a9-6f23-43af-be09-f0b52d3c6310
begin 
	vis_data  = readdlm(joinpath(ethalons_dir , "BaSO4.txt"))
	gold_data_nir = readdlm(joinpath(ethalons_dir , "gold_nir.txt"))
	gold_data_ir = readdlm(joinpath(ethalons_dir , "gold_ir.txt"))
end;

# ╔═╡ 8425182a-b9c4-4e3f-ae06-c20e180b126a
function split_and_mean(m)
	# gets matrix input takes first column as l and evaluates r as an average of other columns
	l = vec(m[:,1])
	r = zeros(eltype(l) , length(l))
	M = size(m , 2) - 1
	@show M
	for i in 1:M 
		_V = @view m[: , i  + 1]
		@. r +=_V
	end
	@. r /=M
	return (l , r)
end

# ╔═╡ 457b939f-f980-4519-a277-2f80f842243e
md"""
### Merging spectra 


"""

# ╔═╡ 9aea683f-dbe8-40d4-86b3-4f19158d7651
md" Bernstein polynomial degree $(@bind visnir_bernstein_N Select(1:1:12 , default = 9))"

# ╔═╡ 63f07a8a-ee1c-4490-96f0-1c01307b3b5a
poly_type = ScaledPolynomials.BernsteinSymPoly{visnir_bernstein_N , Float64}

# ╔═╡ c2d0dba6-ebd7-42f9-b417-1bbda8298299


# ╔═╡ 7c9f849d-73ef-4d9c-b2be-2e706d5652ca
md"""
VIS-NIR mearge range,  μm

$(@bind visnir_merge_range RangeSlider(0.1:1e-3:4 , default = 0.9:1e-3:2.0))

"""

# ╔═╡ 1248e540-255a-4669-b50b-dad0a2a9afe4
md"""
VIS-NIR merge plot x_range,  μm

$(@bind visnir_x_merge_range RangeSlider(0.1:1e-2:5 , default= 1.5:1e-2:1.9))

"""

# ╔═╡ 5b7357bb-5b50-4686-8c9b-a0d3f7a460ad
md"""
NIR-IR merge range,  μm

$(@bind nir_ir_merge_range RangeSlider(2.1:1e-2:6 , default = 2.4:1e-2:3.5))

"""

# ╔═╡ 5a11a9f8-9639-4140-85e1-9969bd5a1c70
md"""
NIR-IR merge plot x_range,  μm

$(@bind nirir_x_merge_range RangeSlider(0.1:1e-2:20 , default= 2:1e-2:4))

"""

# ╔═╡ a7c3dca3-dd67-4d87-afed-591c485f6d7c
md" Bernstein polynomial degree $(@bind nir_ir_bernstein_N Select(1:1:12 , default = 9))"

# ╔═╡ fc3a7561-50e9-4294-b1b7-53f4ddb1975f
@bind save_but Button("Save to file")

# ╔═╡ b628823b-2b59-4556-bb77-696d5fb7a301
p_merged[]

# ╔═╡ 40c50676-1382-4693-8d81-15f8692a14a9
function average_by_x(x, y; sorted=true)
    # Map each x to (sum_of_y, count)
    acc = Dict{eltype(x), Tuple{Float64, Int}}()
    
    for (xi, yi) in zip(x, y)
        s, c = get(acc, xi, (0.0, 0))
        acc[xi] = (s + yi, c + 1)
    end

    # Extract unique x values (optionally sorted)
    x_unique = sorted ? sort!(collect(keys(acc))) : collect(keys(acc))
    
    # Compute means
    y_mean = [acc[xi][1] / acc[xi][2] for xi in x_unique]
    
    return x_unique, y_mean
end

# ╔═╡ 66dde037-402f-4d49-adfb-cde03f2e1c0d
function load_ethalons()
	vis_data  = readdlm(joinpath(ethalons_dir , "BaSO4.txt"))
	gold_data_nir = readdlm(joinpath(ethalons_dir , "gold_nir.txt"))
	gold_data_ir = readdlm(joinpath(ethalons_dir , "gold_ir.txt"))
	
	l_vis , r_vis = split_and_mean(vis_data)
	vis = linear_interpolation(l_vis , r_vis)

	l_gold_nir , r_gold_nir = split_and_mean(gold_data_nir)
	l_gold , r_gold =split_and_mean(gold_data_ir)

	l_ir = vcat(l_gold_nir , l_gold)
	r_ir = vcat(r_gold_nir , r_gold)
	i = sortperm(l_ir)
	
	(l , r)  = average_by_x(l_ir[i], r_ir[i]; sorted=true)
	ir = linear_interpolation(l , r)
	return (vis , ir)
end

# ╔═╡ 0130254d-8d1c-4b52-8ae7-8a2613a61f57
(vis_eth , ir_eth)  = load_ethalons()

# ╔═╡ d4159420-ff69-4a86-aab4-75d99a7dccc4
function parse_sample_number(str)
	reg = r"_(\d+\.\d+)"
	return String(first(collect(first(eachmatch(reg , str)))))
end

# ╔═╡ 46a378b8-cf7c-49cd-8350-28b9171666d9
begin 
	names = (:VIS , :NIR , :IR)
	names_str = [String(f) for f in names]
	data = Dict{String, NamedTuple{names , NTuple{3 ,Matrix{Float64}}}}()
	data_cal = Dict{String, NamedTuple{names , NTuple{3 ,Matrix{Float64}}}}()
	data_names = [ parse_sample_number(f)   for f in readdir(data_dir_ir) if isfile(joinpath(data_dir_ir , f))]
	for n in data_names 
		vals = ntuple(3) do i 
			fld = joinpath(data_dir , String(names[i]))
			cnt = contains(n)
			f_i = first(f for f in readdir(fld) if cnt(f))
			full_file = joinpath(fld , f_i)
			is_v = names[i] == :VIS 
			v = readdlm(full_file , (is_v ? ',' : ';') , Float64 )
			v1 = @view v[: , 1]
			v2 = @view v[: , 2]	
			
			@. v2 /= 100.0
			if is_v
				@. v1 /= 1000.0
			else
				@. v1 = 10000.0 / v1
				_inds = sortperm(v1)
				@. v1 = v1[_inds]
				@. v2 = v2[_inds]
			end
			v_cal = calibrate(v , (is_v ? vis_eth : ir_eth) )
			return (v , v_cal)
		end
		data[n] = Tuple(first(v) for v in vals)
		data_cal[n] = Tuple(last(v) for v in vals)
	end
	
end
	#for (fold , d , mult) in zip(("IR" , "NIR" , "VIS") , ())

# ╔═╡ 3d8ed83a-5019-4584-a7ae-0d2d59116238
md""" 
	$(@bind is_cal Select([:cal , :nocal] , default = :cal)) \


	$(@bind sample2plot PlutoUI.MultiSelect(data_names)) \

	$(@bind types PlutoUI.MultiSelect(names_str, default = names_str)) \

	$(@bind x_range PlutoUI.RangeSlider(0:1e-2:20 , default = 0.1:1e-2:17)) \

	$(@bind y_range PlutoUI.RangeSlider(0:1e-3:1.3 , default = 0.05:1e-3:1.2)) 
	"""	

# ╔═╡ ab35716b-5fd9-4e4f-8ff1-b3d8082e5c09
md" Select sample number $(@bind sample2merge PlutoUI.Select(data_names)) " 

# ╔═╡ 655a7200-5d49-4171-be8e-f92efe6ff7b8
if isassigned(fully_merged) && !isempty(fully_merged[])
	save_but
	result_file = save_file(joinpath(results_dir , "$(sample2merge)_res.csv") , filterlist="csv")
	if length(result_file) > 0 
		out = hcat(fully_merged[]  , 1 .- fully_merged[][:,2])
		writedlm(result_file , out)
		fully_merged[] = Matrix{Float64}(undef  , 0 , 0)
	end
end

# ╔═╡ 55213f0f-f282-47d1-835c-9e738313d6a1
function autoscale_y!(p::Plots.Plot; subplot::Int=1, margin::Real=0.05)
    sp = p.subplots[subplot]
    xmin, xmax = xlims(sp)

    ymin, ymax = Inf, -Inf
    for s in sp.series_list
        x, y = s[:x], s[:y]
        for i in eachindex(x, y)
            xi, yi = x[i], y[i]
            if xmin <= xi <= xmax && isfinite(yi)
                ymin = min(ymin, yi)
                ymax = max(ymax, yi)
            end
        end
    end

    if isfinite(ymin) && isfinite(ymax)
        span = ymax - ymin
        pad = span == 0 ? max(abs(ymin) * margin, 1e-6) : span * margin
        ylims!(sp, (ymin - pad, ymax + pad))
    end
    return p
end

# ╔═╡ 5991d57e-9460-4ccc-933d-929d798e16de
begin 
	pl = plot(;plot_common_args...)
	data_plot = is_cal == :cal ? data_cal : data
	for sm in sample2plot
		for t in types
			m = data_plot[sm][Symbol(t)]
			plot!(pl , m[:,1] , m[:,2] , label = String(sm);plot_common_args...)
		end
	end
	xlims!(pl , extrema(x_range))
	autoscale_y!(pl)
end

# ╔═╡ 9d5092d2-d295-42d3-8c12-cfd81293a569
function merge_data(d1, d2, merge_range::NTuple{2,Float64},     polytype)
    poly_obj = polytype()
    N = ScaledPolynomials.parnumber(poly_obj)
    l1, i1 = @views eachcol(d1)
    l2, i2 = @views eachcol(d2)
    left, right = merge_range

    # out-region points flag 
    f1 = @. (l1 < left)
    f2 = @. (l2 > right)

    # in-region points flag
    in_mr1 = @. (left <= l1 <= right) 
    in_mr2 = @. (left <= l2 <= right)
    
    # out-region views
    or_l1, or_i1 = @views (l1[f1], i1[f1])
    or_l2, or_i2 = @views (l2[f2], i2[f2])

    # in-region views
    inr_l1 , inr_i1 = @views (l1[in_mr1], i1[in_mr1])
    inr_l2 , inr_i2 = @views (l2[in_mr2], i2[in_mr2])
    
    x_all = vcat(inr_l1, inr_l2)
    y_all = vcat(inr_i1, inr_i2)
    
    is_high_degree_poly = N > 4
    
    if is_high_degree_poly
        dummy_poly = ScaledPolynomials.BernsteinSymPoly{4 , Float64}()
        bren_max_loc = ScaledPolynomials.bern_max_locations(dummy_poly)
        ScaledPolynomials.scale_ξ_to_x!(bren_max_loc , left , right)

        
        
        r1 = findfirst(in_mr1):findlast(in_mr1)
        r1_ext = max(firstindex(l1), first(r1) - 1):last(r1)
        left_interp = @views linear_interpolation(l1[r1_ext], i1[r1_ext], extrapolation_bc=Line())
        
        r2 = findfirst(in_mr2):findlast(in_mr2)
        r2_ext = first(r2):min(lastindex(l2), last(r2) + 1)
        right_interp = @views linear_interpolation(l2[r2_ext], i2[r2_ext], extrapolation_bc=Line())
        y_dummy = (left_interp(bren_max_loc[1]) , left_interp(bren_max_loc[2]) ,  right_interp(bren_max_loc[3]) , right_interp(bren_max_loc[4]))
        copyto!(dummy_poly.coeffs , y_dummy)
        for _ in 1:(N - 4)
            dummy_poly = ScaledPolynomials.elevate_degree(dummy_poly)
        end
    end
    sp = ScaledPolynomials.ScaledPolynomial(poly_obj, xmin = left, xmax = right)
    ScaledPolynomials.polyfit!(sp, x_all, y_all)

    if !is_high_degree_poly
        sp.poly.coeffs[1] = last(or_i1)
        sp.poly.coeffs[end] = first(or_i2)
    else
        sp.poly.coeffs[1] = dummy_poly.coeffs[1]
        sp.poly.coeffs[2] = dummy_poly.coeffs[2]
        sp.poly.coeffs[end-1] = dummy_poly.coeffs[end-1]
        sp.poly.coeffs[end] = dummy_poly.coeffs[end]
    end

    merge_l = collect(range(left, right, length = 100))
    merge_vals = sp.(merge_l)

    x_cut = vcat(or_l1, or_l2)
    y_cut = vcat(or_i1, or_i2)

    xm = vcat(or_l1, merge_l, or_l2)
    ym = vcat(or_i1, merge_vals, or_i2)

    return ((x_cut, y_cut), (merge_l, merge_vals), (xm, ym), sp)
end

# ╔═╡ 2f3c7d9d-e590-4d63-943f-26ba3b6111d6
if !isempty(sample2merge) 
	d1 = data_cal[sample2merge][:VIS]
	d2 = data_cal[sample2merge][:NIR]
	mrgd = merge_data(d1, d2, extrema(visnir_merge_range), poly_type )
	vis_nir_merged[] = hcat(mrgd[3][1] , mrgd[3][2])
end;

# ╔═╡ b836d439-8dfc-4354-b09b-73fead9cff9b
begin 
	ppp = plot(;plot_common_args...)
	plot!(ppp , d1[:,1] , d1[:,2] ;plot_common_args... , label = "d1")
	plot!(ppp , d2[:,1] , d2[:,2] ;plot_common_args... , label = "d2")
	plot!(ppp , mrgd[1][1] , mrgd[1][2];plot_common_args... , label = "line")
	plot!(ppp , mrgd[2][1] , mrgd[2][2] ;plot_common_args... , fillrange=0, fillalpha=0.2)
	plot!(ppp , mrgd[3][1] , mrgd[3][2] ;plot_common_args... , label = "merged")
	xlims!(ppp , extrema(visnir_x_merge_range))
	autoscale_y!(ppp)
end

# ╔═╡ 8c2f9219-e895-4b23-89d5-8fd07a003c83
if !isempty(sample2merge) 
	poly_type_ir = ScaledPolynomials.BernsteinSymPoly{nir_ir_bernstein_N , Float64}
	d1_ir = !isassigned(vis_nir_merged) ? data_cal[sample2merge][:NIR] : vis_nir_merged[]
	d2_ir = data_cal[sample2merge][:IR]
	mrgd_ir = merge_data(d1_ir, d2_ir, extrema(nir_ir_merge_range), poly_type_ir )
	fully_merged[] = hcat(mrgd_ir[3][1] , mrgd_ir[3][2])
	refresh_full_plot[] = true
end;

# ╔═╡ c0f01dad-85b1-434b-9537-faabde0eef43
begin 
	ppp_ir = plot(;plot_common_args...)
	plot!(ppp_ir , d1_ir[:,1] , d1_ir[:,2] ;plot_common_args... , label = "d1")
	plot!(ppp_ir , d2_ir[:,1] , d2_ir[:,2] ;plot_common_args... , label = "d2")
	plot!(ppp_ir , mrgd_ir[1][1] , mrgd_ir[1][2];plot_common_args... , label = "line")
	plot!(ppp_ir , mrgd_ir[2][1] , mrgd_ir[2][2] ;plot_common_args... , fillrange=0, fillalpha=0.2)
	plot!(ppp_ir , mrgd_ir[3][1] , mrgd_ir[3][2] ;plot_common_args... , label = "merged")
	xlims!(ppp_ir , extrema(nirir_x_merge_range))
	autoscale_y!(ppp_ir)
end

# ╔═╡ Cell order:
# ╠═06272905-c414-4a00-9373-1e216fc08906
# ╠═9681c7ff-7512-4d80-8eee-6af7e7597ad4
# ╠═af55b1c3-8825-46b1-986b-e2af7e57907a
# ╠═f8f60370-e96f-11f1-9a04-21ec278d93a8
# ╠═a31e3058-eac7-4d6e-9e3e-8376aa9543bf
# ╠═6ba449bc-ff95-4bb2-b023-bd072e56c319
# ╟─e5152a42-16c6-4ee2-a3a9-479b417caba3
# ╟─ee742d64-e0d6-40ed-9733-cec0d59031dd
# ╠═66373b23-b865-4b9f-9042-722bf6b9288f
# ╠═ed3276ba-7c3e-42ea-88c6-4d58fe9c1d3a
# ╠═aec51390-534c-4493-b1f8-e139032acbc9
# ╟─46a378b8-cf7c-49cd-8350-28b9171666d9
# ╟─3cf28b79-09e3-4bed-9b9c-2811fbcdf4d7
# ╟─3d8ed83a-5019-4584-a7ae-0d2d59116238
# ╠═5991d57e-9460-4ccc-933d-929d798e16de
# ╟─ab6c19a9-6f23-43af-be09-f0b52d3c6310
# ╠═0130254d-8d1c-4b52-8ae7-8a2613a61f57
# ╟─8425182a-b9c4-4e3f-ae06-c20e180b126a
# ╟─457b939f-f980-4519-a277-2f80f842243e
# ╟─ab35716b-5fd9-4e4f-8ff1-b3d8082e5c09
# ╟─9aea683f-dbe8-40d4-86b3-4f19158d7651
# ╟─63f07a8a-ee1c-4490-96f0-1c01307b3b5a
# ╟─c2d0dba6-ebd7-42f9-b417-1bbda8298299
# ╟─2f3c7d9d-e590-4d63-943f-26ba3b6111d6
# ╟─7c9f849d-73ef-4d9c-b2be-2e706d5652ca
# ╟─1248e540-255a-4669-b50b-dad0a2a9afe4
# ╟─b836d439-8dfc-4354-b09b-73fead9cff9b
# ╟─5b7357bb-5b50-4686-8c9b-a0d3f7a460ad
# ╟─c0f01dad-85b1-434b-9537-faabde0eef43
# ╟─5a11a9f8-9639-4140-85e1-9969bd5a1c70
# ╟─a7c3dca3-dd67-4d87-afed-591c485f6d7c
# ╟─8c2f9219-e895-4b23-89d5-8fd07a003c83
# ╟─fc3a7561-50e9-4294-b1b7-53f4ddb1975f
# ╠═b628823b-2b59-4556-bb77-696d5fb7a301
# ╟─655a7200-5d49-4171-be8e-f92efe6ff7b8
# ╟─66dde037-402f-4d49-adfb-cde03f2e1c0d
# ╟─40c50676-1382-4693-8d81-15f8692a14a9
# ╟─d4159420-ff69-4a86-aab4-75d99a7dccc4
# ╟─55213f0f-f282-47d1-835c-9e738313d6a1
# ╟─9d5092d2-d295-42d3-8c12-cfd81293a569
