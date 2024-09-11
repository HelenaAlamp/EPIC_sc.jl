# Open and plot .txt file
file1 = "tunex_file.txt"
M1 = readdlm(file1)

file2 = "tuney_file.txt"
M2 = readdlm(file2)


plot(M1[:,1], (M2[:,1]), seriestype=:scatter, markersize=0.5, xlabel = "Δνx", ylabel = "Δνy")



#save in .txt
open("tunex_file.txt", "w") do io
    writedlm(io, Δνx)

open("tuney_file.txt", "w") do io
    writedlm(io, Δνy)
end



#plot from .txt NO WORK
file = open("tunex_file.txt", "r")
data_lines = readline.(file)
close(file)

data = map(Float64, data_lines)

##########################################################

using XLSX
XLSX.openxlsx("data.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "new_sheet")
    sheet["A1"] = "this"
    sheet["A2"] = "is a"
    sheet["A3"] = "new file"
    sheet["A4"] = 100

    # will add a row from "A5" to "E5"
    sheet["A5"] = collect(1:5) # equivalent to `sheet["A5", dim=2] = collect(1:4)`

    # will add a column from "B1" to "B4"
    sheet["B1", dim=1] = collect(1:4)

    # will add a matrix from "A7" to "C9"
    sheet["A7:C9"] = [ 1 2 3 ; 4 5 6 ; 7 8 9 ]

end


data = M[:,3]
csvfile = open("data.csv", "w")
CSV.write(csvfile, data)
close(file)