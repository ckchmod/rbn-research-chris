clear all; clc;
%Metropolis Algorithm
dim_size = 100;
A_disordered = zeros(dim_size,dim_size);
E_disordered_0 = zeros(dim_size,dim_size);
simulation_time = 1000000;
loop_counter = 1;
temp = 2.2:2.2:2.2;
mag = zeros(1,length(temp));
dt = 10000;
E_average_magnetization= zeros(simulation_time/dt + 1,1);
movie_counter = 1;
plotting = true;

%Instantiate matrices: A_disordered and A_disordered
for i = 1:size(A_disordered,1)
    for j = 1:size(A_disordered,2)
        A_disordered(i,j) = randi(2)-1;
        if A_disordered(i,j) == 0
            A_disordered(i,j) = -1;
        end
    end
end

for tt = temp
    temperature = tt;
    % Initial Energy Calculation
    for i = 1:size(A_disordered,1)
        for j = 1:size(A_disordered,2)
            [N,S,W,E] = boundaries(i,j,A_disordered);
            E_disordered_0(i,j) =  -A_disordered(i,j)* (A_disordered(S,j) + A_disordered(N,j) + A_disordered(i,E) + A_disordered(i,W));
        end
    end
    E_average_magnetization = sum(A_disordered, 'all')/(size(A_disordered,1) * size(A_disordered,2));

    if (plotting == true)
        figure(1);
        subplot(1,2,1);
        imagesc(A_disordered), title('Spin States'), colorbar, colormap winter;
        subplot(1,2,2);
        imagesc(E_disordered_0), title('Energy'), colorbar, colormap winter;
        title_time = strcat('t=1 ', ' temperature: ', ' ', string( temperature ));
        suptitle(title_time);
    end

    for t = 1:simulation_time
        mc_i = randi(size(A_disordered,1)); mc_j = randi(size(A_disordered,2));
        [N,S,W,E] = boundaries(mc_i, mc_j, A_disordered);
        E_disordered_new =  -1*(-A_disordered(mc_i,mc_j))* (A_disordered(S,mc_j) + A_disordered(N,mc_j) + A_disordered(mc_i,E) + A_disordered(mc_i,W));
        if (E_disordered_new < E_disordered_0(mc_i,mc_j))
            A_disordered(mc_i,mc_j) = A_disordered(mc_i,mc_j) * -1;
            %E_disordered_0(mc_i,mc_j) = E_disordered_new;
        else
            if (rand() < exp(-(2*E_disordered_new)/temperature)) %(2*E_disordered_new)/temperature)
                A_disordered(mc_i,mc_j) = A_disordered(mc_i,mc_j) * -1;
            end
        end
        if ( mod(t, dt) == 0 && (plotting == true) )
            movie_counter = movie_counter + 1;
            subplot(1,2,1);
            imagesc(A_disordered), title('Spin States'), colorbar, colormap winter;
            subplot(1,2,2);
            imagesc(E_disordered_0), title('Energy'), colorbar, colormap winter;
            title_time = strcat('t= ', string(t), ' temperature: ', ' ', string( temperature ) );
            suptitle(title_time);
            drawnow
        end
        E_average_magnetization = E_average_magnetization + sum(A_disordered, 'all')/(size(A_disordered,1) * size(A_disordered,2));
    end
    E_average_magnetization = E_average_magnetization/(simulation_time+1);
    mag(loop_counter) = E_average_magnetization;
    E_average_magnetization = 0;
    loop_counter = loop_counter+1;
end

figure(2)
plot(temp,mag)
xlabel('Temp'); ylabel('Magnetization')
ylim([-1,1]);
title('Temp vs. Average Magnetization')       