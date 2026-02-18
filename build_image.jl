using PackageCompiler

create_sysimage(
    [:JosephsonCircuits, :Symbolics, :Plots, :DelimitedFiles, :LinearAlgebra, :Distributed];
    sysimage_path="TWPA_sysimg.so",  # 出力ファイル名
    precompile_execution_file="precompile_twpa.jl" # 手順2で作ったファイル
)