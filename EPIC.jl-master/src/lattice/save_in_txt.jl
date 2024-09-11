# function to save in file .txt

function save_to_file(filename, vector::Vector{Float64})
    open(filename, "w") do file
        for elem in vector
            write(file, "$elem\n")
        end
    end
end
