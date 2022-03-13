input = Vector{String}()

open("input.txt") do f
    line = 0  
    while ! eof(f) 
        
       s = readline(f)
       push!(input,s)         
       line += 1
       
    end
   
  end

k = 1
n = size(input)[1] + k

function parseToBinary(inputlist)
    outputmatrix = Vector{Int64}[]
    for i in inputlist
        outputleft = Vector{Int64}()
        outputright = Vector{Int64}()
        for j in i
            binarytuple = pauliBinary(j)
            push!(outputleft,binarytuple[1])
            push!(outputright,binarytuple[2])
        end
        outputrow = append!(outputleft,outputright)
        push!(outputmatrix,outputrow)
    end
    return outputmatrix
end

function pauliBinary(pchar)
    if pchar == 'I'
        return (0,0)
    elseif pchar == 'X'
        return (1,0)
    elseif pchar == 'Z'
        return (0,1)
    elseif pchar == 'Y'
        return (1,1)
    end
end

print(parseToBinary(input))

