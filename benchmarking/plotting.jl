using BenchmarkTools, Plots, LaTeXStrings, Statistics
pgfplotsx()

function reduction_methods()
    fast_ci, standard_ci = [], []
    x_coords = []
    fast_coords, standard_coords = [], []

    for size in [UInt8, UInt16, UInt32, UInt64, UInt128]
        open("benchmarking/reduction_methods/reduce_fast$size.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            fast_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            fast_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end
        open("benchmarking/reduction_methods/reduce_standard$size.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            standard_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            standard_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        p = plot(x_coords, standard_coords, yerror=standard_ci,
            title="$size",
            size = (300,250),
            y_lims = (0,35),
            label="Standard reduce",
            xlabel=L"\log_2 \textrm{field size}",
            ylabel=L"\textrm{time} / \upmu\textrm{s}")

        plot!(p, x_coords, fast_coords, yerror=fast_ci,
            legend= size==UInt8,
            label="Fast reduce")

        savefig("benchmarking/reduction_methods/reduction_methods$size.tex")
    end
end

function shiftandadd_vs_comb()
    shift_ci, comb_ci = [], []
    x_coords = []
    shift_coords, comb_coords = [], []

    for size in [UInt8, UInt16, UInt32, UInt64, UInt128]
        open("benchmarking/shiftandadd_vs_comb/shiftandadd_vs_comb$size.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            comb_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            comb_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
            shift_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            shift_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        p = plot(x_coords, shift_coords, yerror=shift_ci,
            title="$size",
            size = (300,250),
            y_lims = size==UInt8 ? (0,90) : (0, 60),
            label="Shift-and-add method",
            xlabel=L"\log_2 \textrm{field size}",
            ylabel=L"\textrm{time} / \upmu\textrm{s}")

        plot!(p, x_coords, comb_coords, yerror=comb_ci,
            legend= size==UInt8,
            label="Comb method")

        savefig("benchmarking/shiftandadd_vs_comb/shiftandadd_vs_comb$size.tex")
        #savefig("benchmarking/shiftandadd_vs_comb/shiftandadd_vs_comb$size")
    end
end

function window(w, m)
    if m==1 return m/2 end
    return (2^w *w)/2 + (m/w)*(1-0.5^w)
end
function windowsize_fieldmult()
    x_coords = [100, 600]
    p = plot()
    for w in [1,2,4,8]
        y_coords = [window(w,m) for m in x_coords]
        plot!(p, x_coords, y_coords,
        size = (300,200),
        label="w=$w",
        xlabel=L"\log_2 \textrm{field size}",
        ylabel=L"t")
    end
    savefig("benchmarking/windowsize_fieldmult/windowsize_model.tex")
    #savefig("benchmarking/windowsize_fieldmult/windowsize_model")

    x_coords = []
    comb_coords, shift_coords, comb_ci, shift_ci = Dict(), Dict(), Dict(), Dict()
    open("benchmarking/windowsize_fieldmult/fieldmult_windowsizes.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        for w in [1,2,4,8]
            comb_coords[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
            comb_ci[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
            shift_coords[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
            shift_ci[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end
    end

    p = plot(x_coords, shift_coords[1], yerror=shift_ci[1],
    title = "Shift-and-add method",
    y_lims = (0,60),
    size = (300,200),
    legend=false,
    label="w=1",
    xlabel=L"\log_2 \textrm{field size}",
    ylabel=L"\textrm{time} / \upmu\textrm{s}")
    for w in [2,4,8]
        plot!(x_coords, shift_coords[w], yerror=shift_ci[w], label="w=$w")
    end
    savefig("benchmarking/windowsize_fieldmult/windowsize_shiftandadd.tex")
    #savefig("benchmarking/windowsize_fieldmult/windowsize_shiftandadd")

    p = plot(x_coords, comb_coords[1], yerror=comb_ci[1],
    title = "Comb method",
    y_lims = (0,60),
    legend=false,
    size = (300,200),
    label="w=1",
    xlabel=L"\log_2 \textrm{field size}",
    ylabel=L"\textrm{time} / \upmu\textrm{s}")
    for w in [2,4,8]
        plot!(x_coords, comb_coords[w], yerror=comb_ci[w], label="w=$w")
    end
    savefig("benchmarking/windowsize_fieldmult/windowsize_comb.tex")
    #savefig("benchmarking/windowsize_fieldmult/windowsize_comb")
end

function threads(w=1)
    threads_ci, standard_ci = [], []
    x_coords = []
    threads_coords, standard_coords = [], []

    loc = w==1 ? "threads" : "threads_window"

    open("benchmarking/threads/$loc.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        standard_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
        standard_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        threads_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
        threads_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
    end

    p = plot(x_coords, threads_coords, yerror=threads_ci,
        size = (240,200),
        y_lims = (0,30),
        title = w==1 ? "Standard multiplication" : "Multiplication, w=$w",
        label= w==1  ? nothing : "Double-threaded",
        xlabel=L"\log_2 \textrm{field size}",
        ylabel=L"\textrm{time} / \upmu\textrm{s}")

    plot!(p, x_coords, standard_coords, yerror=standard_ci,
        legend= true,
        label= w==1  ? nothing : "Single-threaded")

    savefig("benchmarking/threads/$loc.tex")
    #savefig("benchmarking/$loc/$loc")
end

function windowsize_scalarmult_dicts(loc="mult_standard", ws=[1,2,4,8])
    groups = [SECT163K1, SECT233K1, SECT283K1, SECT409K1, SECT571K1]
    groups = [group(UInt64) for group in groups]
    x_coords = [log(2, group.n) for group in groups]
    y_coords, ci = Dict(), Dict()
    open("benchmarking/$loc/$loc.txt", "r") do io
        eval(Meta.parse(readline(io)))
        for w in ws
            y_coords[w] = [x/1000000 for x in eval(Meta.parse(readline(io)))]
            ci[w] = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        end
    end
    return (x_coords, y_coords, ci)
end
function windowsize_scalarmult(loc="mult_standard", title="Double-and-add Windowsizes", ws=[1,2,4,8], wsplot=[1,2,4,8])
    (x_coords, y_coords, ci) = windowsize_scalarmult_dicts(loc, ws)

    p = plot(x_coords, y_coords[1], yerror=ci[1],
    #title = "$title",
    size = (300,200),
    legend=true,
    label="w=1",
    xlabel=L"\log_2 \textrm{group size}",
    ylabel=L"\textrm{time} / \textrm{ms}")
    for w in wsplot[2:length(wsplot)]
        if w!=1 plot!(x_coords, y_coords[w], yerror=ci[w], label="w=$w") end
    end
    savefig("benchmarking/$loc/$loc.tex")
    savefig("benchmarking/$loc/$loc")
end

function best_scalarmult(w1=4, w2=8, w3=6)
    (x_coords, daa_y, daa_ci) = windowsize_scalarmult_dicts()
    (x_coords, bnaf_y, bnaf_ci) = windowsize_scalarmult_dicts("mult_bnaf")
    (x_coords, wnaf_y, wnaf_ci) = windowsize_scalarmult_dicts("mult_wnaf", [1,2,3,4,5,6,7,8])
    (x_coords, threaded_y, threaded_ci) = windowsize_scalarmult_dicts("mult_threaded", [1])

    p = plot(x_coords, daa_y[w1], yerror=daa_ci[w1],
    size = (300,200),
    legend=true,
    label="Double-and-add, w=4",
    xlabel=L"\log_2 \textrm{group size}",
    ylabel=L"\textrm{time} / \textrm{ms}")

    plot!(x_coords, bnaf_y[w2], yerror=bnaf_ci[w2], label="Binary NAF, w=$w2")
    plot!(x_coords, wnaf_y[w3], yerror=wnaf_ci[w3], label="Width-$w3 NAF")
    plot!(x_coords, threaded_y[1], yerror=threaded_ci[1], label="Threaded binary NAF")

    savefig("benchmarking/scalar_mult/scalar_mult.tex")
    savefig("benchmarking/scalar_mult/scalar_mult")
end

function plot_memo()
    (x_coords, wnaf_y, wnaf_ci) = windowsize_scalarmult_dicts("mult_wnaf", [1,2,3,4,5,6,7,8])
    (x_coords, memo_y, memo_ci) = windowsize_scalarmult_dicts("mult_memo", [1])

    p = plot(x_coords, memo_y[1], yerror=memo_ci[1],
    size = (300,200),
    legend=true,
    label="Memoised multiplication",
    xlabel=L"\log_2 \textrm{group size}",
    ylabel=L"\textrm{time} / \textrm{ms}")

    plot!(x_coords, wnaf_y[6], yerror=wnaf_ci[6], label="Standard multiplication")

    savefig("benchmarking/mult_memo/mult_memo.tex")
    savefig("benchmarking/mult_memo/mult_memo")
end

function double_threads()
    threads_ci, standard_ci = [], []
    x_coords = []
    threads_coords, standard_coords = [], []

    open("benchmarking/threads_double/threads_double.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        standard_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
        standard_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        threads_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
        threads_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
    end

    p = plot(x_coords, threads_coords, yerror=threads_ci,
        size = (300,200),
        label= "Multithreaded",
        xlabel=L"\log_2 \textrm{group order}",
        ylabel=L"\textrm{time} / \upmu\textrm{s}")

    plot!(p, x_coords, standard_coords, yerror=standard_ci,
        legend= true,
        label= "Single threaded")

    savefig("benchmarking/threads_double/threads_double.tex")
    #savefig("benchmarking/threads_double/threads_double")
end

function mult_mont()
    threads_ci, standard_ci = [], []
    x_coords = []
    threads_coords, standard_coords = [], []

    open("benchmarking/mult_mont/mult_mont.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        standard_coords = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        standard_ci = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        threads_coords = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        threads_ci = [x/1000000 for x in eval(Meta.parse(readline(io)))]
    end

    (x_coords, wnaf_y, wnaf_ci) = windowsize_scalarmult_dicts("mult_threaded", [1])

    p = plot(x_coords, threads_coords, yerror=threads_ci,
        size = (300,200),
        label= "Affine Montgomery powering ladder",
        xlabel=L"\log_2 \textrm{group size}",
        ylabel=L"\textrm{time} / \textrm{ms}")

    plot!(p, x_coords, standard_coords, yerror=standard_ci,
        legend= true,
        label= "General Montgomery powering ladder")

    plot!(p, x_coords, wnaf_y[1], yerror=wnaf_ci[1],
        legend= true,
        label= "Multithreaded binary NAF method")

    savefig("benchmarking/mult_mont/mult_mont.tex")
    savefig("benchmarking/mult_mont/mult_mont")
end

function plot_packages()

    for (op, title) in [("mult", "Multiplication"), ("sq", "Squaring"), ("inv", "Inversion")]
        becc_ci, other_ci = [], []
        x_coords = []
        becc_coords, other_coords = [], []

        open("benchmarking/becc/becc-$op.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            becc_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            becc_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        open("benchmarking/gf/gf-$op.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            other_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            other_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        #r = 1:length(x_coords)
        p = plot(x_coords, other_coords, yerror=other_ci,
            size = (300,200),
            title = title,
            label= "GaloisFields",
            xlabel=L"\log_2 \textrm{field size}",
            ylabel=L"\textrm{time} / \upmu\textrm{s}")

        open("benchmarking/nemo/nemo-$op.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            other_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            other_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        plot!(p, x_coords, other_coords, yerror=other_ci,
            legend= op=="sq",
            label= "Nemo")

        plot!(p, x_coords, becc_coords, yerror=becc_ci,
            legend= op=="sq",
            label= "BinaryECC")

        savefig("benchmarking/bfield/bfield-$op.tex")
        #savefig("benchmarking/$package/$package-$op")
    end
end

function openssl()
    x_coords, y_coords, ci = [], Dict(), Dict()
    open("benchmarking/openssl/openssl.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        y_coords = eval(Meta.parse(readline(io)))
        ci = eval(Meta.parse(readline(io)))
    end

    for type in ["becc", "openssl"]
        y_coords[type] = [t/1000000 for t in y_coords[type]]
        ci[type] = [t/1000000 for t in ci[type]]
    end

    p = plot(x_coords, y_coords["openssl"], yerror=ci["openssl"],
        size = (300,200),
        label= "OpenSSL",
        xlabel=L"\log_2 \textrm{group size}",
        ylabel=L"\textrm{time} / \textrm{ms}")

    threads_coords, threads_ci = [], []
    open("benchmarking/mult_mont/mult_mont.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        standard_coords = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        standard_ci = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        threads_coords = [x/1000000 for x in eval(Meta.parse(readline(io)))]
        threads_ci = [x/1000000 for x in eval(Meta.parse(readline(io)))]
    end

    plot!(p, x_coords,threads_coords, yerror=threads_ci,
        label= "BinaryECC - Montgomery")

    plot!(x_coords, y_coords["becc"], yerror=ci["becc"],
        label= "BinaryECC - fast", legend=true)

    savefig("benchmarking/openssl/openssl.tex")
    savefig("benchmarking/openssl/openssl")
end

function plot_timing()
    x_coords = []
    y_coords = Dict("mont1"=>[], "mont0"=>[], "std1"=>[], "std0"=>[])
    ci = Dict("mont1"=>[], "mont0"=>[], "std1"=>[], "std0"=>[])
    open("benchmarking/timing/timing.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        y_coords = eval(Meta.parse(readline(io)))
        ci = eval(Meta.parse(readline(io)))
    end

    for key in keys(y_coords)
        y_coords[key] = [t/1000000 for t in y_coords[key]]
        ci[key] = [t/1000000 for t in ci[key]]
    end


    p = plot(x_coords, y_coords["mont1"] .- y_coords["mont0"],
        yerror=ci["mont1"],
        size = (300,200),
        label= "mont diff",
        xlabel=L"\log_2 \textrm{group size}",
        ylabel=L"\textrm{time} / \textrm{ms}")

    plot!(x_coords, [0 for x in x_coords], yerror=ci["mont0"])

    savefig("benchmarking/timing/mont.tex")
    savefig("benchmarking/timing/mont")

    p = plot(x_coords, y_coords["std1"] .- y_coords["std0"],
        yerror=ci["std1"],
        size = (300,200),
        label= "std diff",
        xlabel=L"\log_2 \textrm{group size}",
        ylabel=L"\textrm{time} / \textrm{ms}")

    plot!(x_coords, [0 for x in x_coords], yerror=ci["std0"])

    savefig("benchmarking/timing/std.tex")
    savefig("benchmarking/timing/std")
end
