using BenchmarkTools, Statistics

#95% confidence interval for gaussian with unkown mu
function gaussian_ci(sd, n)
    return 1.96*(sd / sqrt(n))
end

#reduction method type for each word size
function collect_reduce(type="standard")
    x_coords = [113, 131, 163, 193, 233, 239, 283, 409, 571]
    fields = [BFieldPoint113, BFieldPoint131, BFieldPoint163, BFieldPoint193,
        BFieldPoint233, BFieldPoint239, BFieldPoint283, BFieldPoint409, BFieldPoint571]
    sizes = [UInt8,UInt16,UInt32,UInt64,UInt128]

    for size in sizes
        println("wordsize: $size")
        y_coords, ci = [], []
        for i in 1:length(fields)
            field = fields[i]{size}
            println("field: $field")
            L = ceil(Int, 2*x_coords[i]/bitsize(size))
            x = "@benchmark reduce(\$($field(random(StaticUInt{$L,$size}))))"
            x = Meta.parse(x)
            b = eval(x)
            append!(y_coords, [mean(b.times)])
            append!(ci, [gaussian_ci(std(b.times), length(b.times))])
        end

        open("benchmarking/reduction_methods/reduce_$type$size.txt", "w") do io
            write(io, "$x_coords\n")
            write(io, "$y_coords\n")
            write(io, "$ci\n")
        end
    end
end

#shift-and-add vs comb for each wordsize
function collect_shiftandadd_vs_comb()
    x_coords = [113, 131, 163, 193, 233, 239, 283, 409, 571]
    fields = [BFieldPoint113, BFieldPoint131, BFieldPoint163, BFieldPoint193,
        BFieldPoint233, BFieldPoint239, BFieldPoint283, BFieldPoint409, BFieldPoint571]
    sizes = [UInt8,UInt16,UInt32,UInt64,UInt128]

    for size in sizes
        println("wordsize: $size")
        comb_coords, shift_coords, comb_ci, shift_ci = [], [], [], []
        for field in fields
            field = field{UInt64}
            println("field: $field")
            x = "@benchmark mult_comb_rtl(\$(random($field)), \$(random($field)))"
            x = Meta.parse(x)
            b = eval(x)
            append!(comb_coords, [mean(b.times)])
            append!(comb_ci, [gaussian_ci(std(b.times), length(b.times))])
            x = "@benchmark mult_shiftandadd(\$(random($field)), \$(random($field)))"
            x = Meta.parse(x)
            b = eval(x)
            append!(shift_coords, [mean(b.times)])
            append!(shift_ci, [gaussian_ci(std(b.times), length(b.times))])
        end

        open("benchmarking/shiftandadd_vs_comb/shiftandadd_vs_comb$size.txt", "w") do io
            write(io, "$x_coords\n")
            write(io, "$comb_coords\n")
            write(io, "$comb_ci\n")
            write(io, "$shift_coords\n")
            write(io, "$shift_ci\n")
        end
    end
end

function collect_windowsize_fieldmult()
    x_coords = [113, 131, 163, 193, 233, 239, 283, 409, 571]
    fields = [BFieldPoint113, BFieldPoint131, BFieldPoint163, BFieldPoint193,
        BFieldPoint233, BFieldPoint239, BFieldPoint283, BFieldPoint409, BFieldPoint571]

    comb_coords, shift_coords, comb_ci, shift_ci = Dict(), Dict(), Dict(), Dict()
    for w in [1,2,4,8]
        println("windowsize: $w")
        comb_coords[w], shift_coords[w], comb_ci[w], shift_ci[w] = [],[],[],[]
        for field in fields
            field = field{UInt64}
            println("field: $field")
            x = "@benchmark mult_comb_ltr(\$(random($field)), \$(random($field)), $w)"
            x = Meta.parse(x)
            b = eval(x)
            append!(comb_coords[w], [mean(b.times)])
            append!(comb_ci[w], [gaussian_ci(std(b.times), length(b.times))])
            x = "@benchmark mult_shiftandadd(\$(random($field)), \$(random($field)), $w)"
            x = Meta.parse(x)
            b = eval(x)
            append!(shift_coords[w], [mean(b.times)])
            append!(shift_ci[w], [gaussian_ci(std(b.times), length(b.times))])
        end
    end
    open("benchmarking/windowsize_fieldmult/fieldmult_windowsizes.txt", "w") do io
        write(io, "$x_coords\n")
        for w in [1,2,4,8]
            write(io, "$(comb_coords[w])\n")
            write(io, "$(comb_ci[w])\n")
            write(io, "$(shift_coords[w])\n")
            write(io, "$(shift_ci[w])\n")
        end
    end
end

function collect_threads(w=1)
    x_coords = [113, 131, 163, 193, 233, 239, 283, 409, 571]
    fields = [BFieldPoint113, BFieldPoint131, BFieldPoint163, BFieldPoint193,
        BFieldPoint233, BFieldPoint239, BFieldPoint283, BFieldPoint409, BFieldPoint571]
    fields = [f{UInt} for f in fields]

    standard_coords, threaded_coords, standard_ci, threaded_ci = [],[],[],[]
    for field in fields
        println("field: $field")
        if w==1
            x = "@benchmark mult_shiftandadd(\$(random($field)), \$(random($field)))"
        else
            x = "@benchmark mult_shiftandadd_window(\$(random($field)), \$(random($field)), $w)"
        end
        x = Meta.parse(x)
        b = eval(x)
        append!(standard_coords, [mean(b.times)])
        append!(standard_ci, [gaussian_ci(std(b.times), length(b.times))])
        if w==1
            x = "@benchmark mult_threaded(\$(random($field)), \$(random($field)))"
        else
            x = "@benchmark mult_threaded_window(\$(random($field)), \$(random($field)), $w)"
        end
        x = Meta.parse(x)
        b = eval(x)
        append!(threaded_coords, [mean(b.times)])
        append!(threaded_ci, [gaussian_ci(std(b.times), length(b.times))])
    end

    loc = w==1 ? "threads" : "threads_window"
    open("benchmarking/$loc/$loc.txt", "w") do io
        write(io, "$x_coords\n")
        write(io, "$standard_coords\n")
        write(io, "$standard_ci\n")
        write(io, "$threaded_coords\n")
        write(io, "$threaded_ci\n")
    end
end

function collect_threads_window(w=4)
    x_coords = [113, 131, 163, 193, 233, 239, 283, 409, 571]
    fields = [BFieldPoint113, BFieldPoint131, BFieldPoint163, BFieldPoint193,
        BFieldPoint233, BFieldPoint239, BFieldPoint283, BFieldPoint409, BFieldPoint571]

    standard_coords, threaded_coords, standard_ci, threaded_ci = [],[],[],[]
    for field in fields
        field = field{UInt64}
        println("field: $field")
        x = "@benchmark mult_shiftandadd_window(\$(random($field)), \$(random($field)), $w)"
        x = Meta.parse(x)
        b = eval(x)
        append!(standard_coords, [mean(b.times)])
        append!(standard_ci, [gaussian_ci(std(b.times), length(b.times))])
        x = "@benchmark mult_threaded_window(\$(random($field)), \$(random($field)), $w)"
        x = Meta.parse(x)
        b = eval(x)
        append!(threaded_coords, [mean(b.times)])
        append!(threaded_ci, [gaussian_ci(std(b.times), length(b.times))])
    end

    open("benchmarking/threads_window/threads_window.txt", "w") do io
        write(io, "$x_coords\n")
        write(io, "$standard_coords\n")
        write(io, "$standard_ci\n")
        write(io, "$threaded_coords\n")
        write(io, "$threaded_ci\n")
    end
end

function collect_mult_standard()
    x_coords = [163, 233, 283, 409, 571]
    groups = [SECT163K1, SECT233K1, SECT283K1, SECT409K1, SECT571K1]
    y_coords, ci = Dict(), Dict()

    for w in [1,2,4,8]
        println("w: $w")
        y_coords[w], ci[w] = [], []

        for group in groups
            group = group(UInt64)
            println("group: $group")
            G = group.G
            n = group.n

            x = "@benchmark mult_window(\$($G * rand(1:$n)), \$(rand(1:$n)), $w)"
            x = Meta.parse(x)
            b = eval(x)
            append!(y_coords[w], [mean(b.times)])
            append!(ci[w], [gaussian_ci(std(b.times), length(b.times))])
        end
    end

    open("benchmarking/mult_standard/mult_standard.txt", "w") do io
        write(io, "$x_coords\n")
        for w in [1,2,4,8]
            write(io, "$(y_coords[w])\n")
            write(io, "$(ci[w])\n")
        end
    end
end

function collect_mult_bnaf()
    x_coords = [163, 233, 283, 409, 571]
    groups = [SECT163K1, SECT233K1, SECT283K1, SECT409K1, SECT571K1]
    y_coords, ci = Dict(), Dict()

    for w in [1,2,4,8]
        println("w: $w")
        y_coords[w], ci[w] = [], []

        for group in groups
            group = group(UInt64)
            println("group: $group")
            G = group.G
            n = group.n

            x = "@benchmark mult_bnaf_window(\$($G * rand(1:$n)), \$(rand(1:$n)), $w)"
            x = Meta.parse(x)
            b = eval(x)
            append!(y_coords[w], [mean(b.times)])
            append!(ci[w], [gaussian_ci(std(b.times), length(b.times))])
        end
    end

    open("benchmarking/mult_bnaf/mult_bnaf.txt", "w") do io
        write(io, "$x_coords\n")
        for w in [1,2,4,8]
            write(io, "$(y_coords[w])\n")
            write(io, "$(ci[w])\n")
        end
    end
end

function collect_mult_wnaf()
    groups = [SECT163K1, SECT233K1, SECT283K1, SECT409K1, SECT571K1]
    groups = [group(UInt64) for group in groups]
    x_coords = [log(2, group.n) for group in groups]
    y_coords, ci = Dict(), Dict()

    for w in [1,2,3,4,5,6,7,8]
        println("w: $w")
        y_coords[w], ci[w] = [], []

        for group in groups
            println("group: $group")
            G = group.G
            n = group.n

            x = "@benchmark mult_wnaf(\$($G * rand(1:$n)), \$(rand(1:$n)), $w)"
            x = Meta.parse(x)
            b = eval(x)
            append!(y_coords[w], [mean(b.times)])
            append!(ci[w], [gaussian_ci(std(b.times), length(b.times))])
        end
    end

    open("benchmarking/mult_wnaf/mult_wnaf.txt", "w") do io
        write(io, "$x_coords\n")
        for w in [1,2,3,4,5,6,7,8]
            write(io, "$(y_coords[w])\n")
            write(io, "$(ci[w])\n")
        end
    end
end

function collect_double_threaded()
    groups = [SECT163K1, SECT233K1, SECT283K1, SECT409K1, SECT571K1]
    groups = [group(UInt64) for group in groups]
    x_coords = [log(2, group.n) for group in groups]

    standard_coords, threaded_coords, standard_ci, threaded_ci = [],[],[],[]
    for group in groups
        println("group: $group\n")
        G = group.G
        n = group.n
        x = "@benchmark double( \$($G * rand(1:$n)) )"
        x = Meta.parse(x)
        b = eval(x)
        append!(standard_coords, [mean(b.times)])
        append!(standard_ci, [gaussian_ci(std(b.times), length(b.times))])
        x = "@benchmark double_threaded( \$($G * rand(1:$n)) )"
        x = Meta.parse(x)
        b = eval(x)
        append!(threaded_coords, [mean(b.times)])
        append!(threaded_ci, [gaussian_ci(std(b.times), length(b.times))])
    end

    open("benchmarking/threads_double/threads_double.txt", "w") do io
        write(io, "$x_coords\n")
        write(io, "$standard_coords\n")
        write(io, "$standard_ci\n")
        write(io, "$threaded_coords\n")
        write(io, "$threaded_ci\n")
    end
end
