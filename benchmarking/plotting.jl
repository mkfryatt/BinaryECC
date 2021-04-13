using BenchmarkTools, Plots, LaTeXStrings, Statistics
pgfplotsx()

function reduction_methods()
    fast_ci, standard_ci = [], []
    x_coords = []
    fast_coords, standard_coords = [], []

    for size in [8, 16, 32, 64, 128]
        open("benchmarking\\reduction_methods\\reduce_fast$size.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            fast_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            fast_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end
        open("benchmarking\\reduction_methods\\reduce_standard$size.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            standard_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            standard_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        p = plot(x_coords, fast_coords, yerror=fast_ci,
            title="UInt$size",
            size = (300,200),
            label="Fast reduce",
            xlabel=L"\log_2 \textrm{field size}",
            ylabel=L"\textrm{time} / \mu s")

        plot!(p, x_coords, standard_coords, yerror=standard_ci,
            legend= size==8,
            label="Standard reduce")

        savefig("benchmarking\\reduction_methods$size.tex")
    end
end

function shiftandadd_vs_comb()
    shift_ci, comb_ci = [], []
    x_coords = []
    shift_coords, comb_coords = [], []

    for size in [8, 16, 32, 64, 128]
        open("benchmarking\\shiftandadd_vs_comb\\shiftandadd_vs_comb$size.txt", "r") do io
            x_coords = eval(Meta.parse(readline(io)))
            shift_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            comb_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
            shift_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
            comb_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end

        p = plot(x_coords, shift_coords, yerror=shift_ci,
            title="UInt$size",
            size = (300,200),
            label="Shift-and-add method",
            xlabel=L"\log_2 \textrm{field size}",
            ylabel=L"\textrm{time} / \mu s")

        plot!(p, x_coords, comb_coords, yerror=comb_ci,
            legend= size==8,
            label="Comb method")

        savefig("benchmarking\\shiftandadd_vs_comb$size.tex")
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
    savefig("benchmarking\\windowsize_model.tex")

        y_coords = Dict()
        ci = Dict()
    open("benchmarking\\windowsize_fieldmult\\right_to_left_mult-windowsizes.txt") do io
        x_coords = eval(Meta.parse(readline(io)))
        for w in [1,2,4,8]
            y_coords[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
            ci[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end
    end

    p = plot(x_coords, y_coords[1], yerror=ci[1],
    title = "Shift-and-add method",
    size = (300,200),
    legend=false,
    label="w=1",
    xlabel=L"\log_2 \textrm{field size}",
    ylabel=L"\textrm{time} / \mu s")
    for w in [2,4,8]
        plot!(x_coords, y_coords[w], yerror=ci[w], label="w=$w")
    end
    savefig("benchmarking\\windowsize_shiftandadd.tex")

    y_coords = Dict()
    ci = Dict()
    open("benchmarking\\windowsize_fieldmult\\left_to_right_comb_mult-windowsizes.txt") do io
        x_coords = eval(Meta.parse(readline(io)))
        y_coords = Dict()
        ci = Dict()
        for w in [1,2,4,8]
            y_coords[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
            ci[w] = [x/1000 for x in eval(Meta.parse(readline(io)))]
        end
    end

    p = plot(x_coords, y_coords[1], yerror=ci[1],
    title = "Comb method",
    legend=false,
    size = (300,200),
    label="w=1",
    xlabel=L"\log_2 \textrm{field size}",
    ylabel=L"\textrm{time} / \mu s")
    for w in [2,4,8]
        plot!(x_coords, y_coords[w], yerror=ci[w], label="w=$w")
    end
    savefig("benchmarking\\windowsize_comb.tex")


    return p
end

function threads()
    threads_ci, standard_ci = [], []
    x_coords = []
    threads_coords, standard_coords = [], []

    open("benchmarking\\threads\\threads.txt", "r") do io
        x_coords = eval(Meta.parse(readline(io)))
        threads_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
        standard_coords = [x/1000 for x in eval(Meta.parse(readline(io)))]
        threads_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
        standard_ci = [x/1000 for x in eval(Meta.parse(readline(io)))]
    end

    p = plot(x_coords, threads_coords, yerror=threads_ci,
        size = (300,200),
        label="Multithreaded",
        xlabel=L"\log_2 \textrm{field size}",
        ylabel=L"\textrm{time} / \mu s")

    plot!(p, x_coords, standard_coords, yerror=standard_ci,
        legend= true,
        label="Single threaded")

    savefig("benchmarking\\threads.tex")
    savefig("benchmarking\\threads")
end
