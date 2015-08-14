% Example: the running time of the spherical needlet transform

B = 2;

rep = 10;
time = zeros(8, rep);

for power = 2:9
    l_max = 2^power;
    alm = zeros(l_max+1, 2*l_max+1);
    for r = 1:rep
        for l = 1:l_max
            for m = 0:l
                alm(l+1, m+l_max+1) = randn+1i*randn;
            end
        end
        tic
        beta = spneedlet_tran( alm, l_max, B );
        time(power-1, r) = toc;
    end
end

save('run_time.mat', 'time')

boxplot(time', 'labels', {'4', '8', '16', '32', '64', '128', '256', '512'})
med = median(time, 2);
hold on
plot(1:8, med, '--o', 'LineWidth', 1.5)
xlabel('l_{max}', 'FontSize', 15)
ylabel('Time (s)', 'FontSize', 15)
axis tight
set(gca, 'YScale', 'log') 