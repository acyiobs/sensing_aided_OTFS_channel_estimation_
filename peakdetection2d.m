function peak = peakdetection2d(map, train_cell, guard_cell, threshold_factor)

[Nr, Nc] = size(map);

Tr = train_cell(1);
Tc = train_cell(2);

Gr = guard_cell(1);
Gc = guard_cell(2);

peak = zeros(size(map));

for i = Tr+Gr+1:(Nr-Tr-Gr)
    for j = Tc+Gc+1:(Nc-Tc-Gc)
        noise_level = zeros(1,1);
        for k = (i-Tr-Gr) : (i+Tr+Gr)
            for h = (j-Tc-Gc) : (j+Gc+Tc)
                if(abs(k-i)>Gr||abs(h-j)>Gc)
                    noise_level = noise_level + map(k,h);
                end
            end
        end
        length = 2*(Tr+Gr)+1;
        width =  2*(Tc+Gc)+1;
        threshold = noise_level/(length*width-(2*Gr+1)*(2*Gc+1))*threshold_factor;
        peak(i,j) = (map(i,j)>threshold) * (map(i,j)/threshold);
    end
end

end