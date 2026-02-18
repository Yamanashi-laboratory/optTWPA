using JosephsonCircuits
using Symbolics
using Plots
using DelimitedFiles
using LinearAlgebra

# 行列演算のスレッド競合を防ぎ、Juliaのスレッド並列に計算リソースを集中させる
BLAS.set_num_threads(1)

"""
メインの計算処理
"""
function main()
    # 1. 回路変数の定義
    @variables Rleft Rright Cg Lj Cj Cc Cr Lr Cn
    circuit = Tuple{String,String,String,Num}[]

    # 回路構成の組み立て (Nj=2400の大規模回路)
    push!(circuit, ("P$(1)_$(0)", "1", "0", 1))
    push!(circuit, ("R$(1)_$(0)", "1", "0", Rleft))
    push!(circuit, ("C$(1)_$(0)", "1", "0", Cg))
    push!(circuit, ("Lj$(1)_$(2)", "1", "2", Lj))
    push!(circuit, ("C$(1)_$(2)", "1", "2", Cj))

    Nj = 2400
    pmrpitch = 3
    j = 2
    for i in 2:Nj-1
        if mod(i, pmrpitch) == pmrpitch ÷ 3
            push!(circuit, ("C$(j)_$(0)", "$(j)", "0", Cn))
            push!(circuit, ("Lj$(j)_$(j+2)", "$(j)", "$(j+2)", Lj))
            push!(circuit, ("C$(j)_$(j+2)", "$(j)", "$(j+2)", Cj))

            push!(circuit, ("C$(j)_$(j+1)", "$(j)", "$(j+1)", Cc))
            push!(circuit, ("C$(j+1)_$(0)", "$(j+1)", "0", Cr))
            push!(circuit, ("L$(j+1)_$(0)", "$(j+1)", "0", Lr))
            j += 1
        else
            push!(circuit, ("C$(j)_$(0)", "$(j)", "0", Cg))
            push!(circuit, ("Lj$(j)_$(j+1)", "$(j)", "$(j+1)", Lj))
            push!(circuit, ("C$(j)_$(j+1)", "$(j)", "$(j+1)", Cj))
        end
        j += 1
    end

    push!(circuit, ("R$(j)_$(0)", "$(j)", "0", Rright))
    push!(circuit, ("P$(j)_$(0)", "$(j)", "0", 2))

    circuitdefs = Dict(
        Lj => IctoLj(9e-06),
        Cg => 11.0647e-15,
        Cc => 14.1185e-15,
        Cn => 9.44222e-15,
        Cr => 2.5e-12,
        Lr => 130e-12,
        Cj => 369e-15,
        Rleft => 50.0,
        Rright => 50.0,
    )
    ws  = 2π * (1.0:0.01:16.0) * 1e9    #ここで周波数(横軸を変更): (開始値:ステップ幅:終了値)
    wp  = (2π*8.6938e9,)
    Ip  = 5.5168e-06
    sources = [(mode=(1,), port=1, current=Ip)]
    Npumpharmonics       = (20,)
    Nmodulationharmonics = (10,)

    # 2. 並列計算の実行
    nt = Threads.nthreads()
    #println("利用可能なスレッド数: $nt")
    if nt == 1
        @warn "現在1スレッドで動作しています。高速化するには環境変数 JULIA_NUM_THREADS を設定してください。"
    end

    len = length(ws)
    # 各スレッドの結果を格納する配列を準備
    # (通信コストが発生しないため、プロセス並列より高速)
    s21_results = Vector{Vector{ComplexF64}}(undef, nt)
    
    #println("計算開始（hbsolveをマルチスレッドで実行中）...")
    @time Threads.@threads for i in 1:nt
        # 各スレッドの担当範囲を計算
        idx_start = ((i-1) * len ÷ nt) + 1
        idx_end   = (i * len ÷ nt)
        ws_part   = ws[idx_start:idx_end]
        
        # 各スレッドで独立してhbsolveを実行
        res = hbsolve(ws_part, wp, sources, Nmodulationharmonics, 
                      Npumpharmonics, circuit, circuitdefs)
        
        # 必要なS21のみを抽出
        s21_results[i] = res.linearized.S((0,), 2, (0,), 1, :)
    end

    # 3. 結果の結合と保存
    s21_combined = vcat(s21_results...)
    freq_vec_ghz = collect(ws ./ (2π * 1e9))
    gain_vec_db = 10 .* log10.(abs2.(s21_combined))

    # オリジナルの plot_gain 関数が読み込める形式で保存
    open("freq_gain_sim.csv", "w") do io
        println(io, "frequency_GHz,gain_dB")
        for (f, g) in zip(freq_vec_ghz, gain_vec_db)
            println(io, f, ",", g)
        end
    end
    #println("計算完了。データを 'freq_gain_sim.csv' に保存しました。")
end

"""
オリジナルのプロット関数 (そのまま維持)
"""
function plot_gain()
    if !isfile("freq_gain_sim.csv")
        error("CSVファイルが見つかりません。先に main() を実行してください。")
    end

    # CSV読み込み (ヘッダ行を１行スキップ)
    data = readdlm("freq_gain_sim.csv", ',', skipstart=1)
    x = vec(data[:, 1])   # Frequency [GHz]
    y = vec(data[:, 2])   # Gain [dB]

    # プロット
    plt = plot(
        x, y;
        xlabel     = "Frequency [GHz]",
        ylabel     = "Gain [dB]",
        legend     = false,
        seriestype = :line,
        marker     = :none,
        line       = (:solid, 2),
        ylims      = (-1, maximum(y)*1.1),
    )

    display(plt)
    savefig(plt, "gain_.png")
    #println("グラフを 'gain_.png' に保存しました。")
end

# スクリプト実行
if abspath(PROGRAM_FILE) == @__FILE__
    main()
    plot_gain()
end