%% ================================================================
%  CLIMATOLOGIA MENSAL (ECMWF IFS – 31 ANOS)
%  Autor: Augusto Pereira
%  Última atualização: 2025-11-05
%  ================================================================

clc; clear; close all;

%% ==== 1. Leitura das variáveis comuns (coordenadas e tempo) ====
tpfile  = 'tp.nc';                              % precipitação acumulada diária (avgad)
metfile = 'data_stream-moda_stepType-avgua.nc'; % variáveis meteorológicas mensais
oceanfile = 'height_waves_swell.nc';

lon = ncread(tpfile, 'longitude');   % [°E]
lat = ncread(tpfile, 'latitude');    % [°N]
time = ncread(tpfile, 'valid_time'); % [s desde 1970-01-01]

LON = ncread(oceanfile,'longitude');
LAT = ncread(oceanfile,'latitude');

ref_time = datetime(1970,1,1,0,0,0);
time = ref_time + seconds(time);
[yr, mo, ~] = ymd(time);

anos  = unique(yr);
nAnos = numel(anos);
nMes  = 12;
nLon  = length(lon);
nLat  = length(lat);

NLON = length(LON);
NLAT = length(LAT);

%% ==== 2. Leitura das variáveis ====
tp  = ncread(tpfile, 'tp');   % [m] média de acumulados diários
t2m = ncread(metfile, 't2m'); % [K]
u10 = ncread(metfile, 'u10'); % [m/s]
v10 = ncread(metfile, 'v10'); % [m/s]
sp  = ncread(metfile, 'sp');  % [Pa]

swh = ncread(oceanfile, 'swh'); % [m]


%% ==== 3. Converter precipitação para total mensal (rotina manual) ====
% ECMWF: "tp" (avgad) = média do acumulado diário → multiplica por nº de dias do mês

ndays = zeros(size(mo));
for k = 1:numel(mo)
    m = mo(k);
    y = yr(k);

    if any(m == [1,3,5,7,8,10,12])
        ndays(k) = 31;
    elseif any(m == [4,6,9,11])
        ndays(k) = 30;
    elseif m == 2
        % Ano bissexto → fevereiro = 29 dias
        if mod(y,4)==0 && (mod(y,100)~=0 || mod(y,400)==0)
            ndays(k) = 29;
        else
            ndays(k) = 28;
        end
    end
end

tp_total = tp;  % inicializa matriz
for k = 1:length(time)
    tp_total(:,:,k) = tp(:,:,k) * ndays(k);  % total mensal [m/mês]
end

%% ==== 4. Pré-alocação das climatologias ====
clim_tp  = nan(nLon, nLat, nMes);
clim_t2m = nan(nLon, nLat, nMes);
clim_u10 = nan(nLon, nLat, nMes);
clim_v10 = nan(nLon, nLat, nMes);
clim_sp  = nan(nLon, nLat, nMes);

clima_swh  = nan(NLON, NLAT, nMes); 
%% ==== 5. Cálculo explícito da climatologia mensal ====
for m = 1:nMes
    idx_m = (mo == m);

    soma_tp  = zeros(nLon, nLat);
    soma_t2m = zeros(nLon, nLat);
    soma_u10 = zeros(nLon, nLat);
    soma_v10 = zeros(nLon, nLat);
    soma_sp  = zeros(nLon, nLat);
    soma_swh = zeros(nLon, nLat);

    for y = 1:nAnos
        ano_atual = anos(y);
        idx = (mo == m & yr == ano_atual);
        if isempty(find(idx,1)), continue; end

        soma_tp  = soma_tp  + mean(tp_total(:,:,idx), 3, 'omitnan');
        soma_t2m = soma_t2m + mean(t2m(:,:,idx), 3, 'omitnan');
        soma_u10 = soma_u10 + mean(u10(:,:,idx), 3, 'omitnan');
        soma_v10 = soma_v10 + mean(v10(:,:,idx), 3, 'omitnan');
        soma_sp  = soma_sp  + mean(sp(:,:,idx),  3, 'omitnan');

        soma_swh = soma_swh + mean(swh(:,:,idx),  3, 'omitnan');
    end

    clim_tp(:,:,m)  = soma_tp  / nAnos;  % [m/mês]
    clim_t2m(:,:,m) = soma_t2m / nAnos;  % [K]
    clim_u10(:,:,m) = soma_u10 / nAnos;  % [m/s]
    clim_v10(:,:,m) = soma_v10 / nAnos;  % [m/s]
    clim_sp(:,:,m)  = soma_sp  / nAnos;  % [Pa]
    
    clima_swh(:,:,m)  = soma_swh  / nAnos;  % [m]

    fprintf(' Climatologia do mês %02d concluída (%d anos)\n', m, nAnos);
end

%% ==== 6. Conversões opcionais (visualização) ====
 t2m_C  = clim_t2m - 273.15;  % [°C]
 tp_mm  = clim_tp * 1000;     % [mm/mês]
 sp_hPa = clim_sp / 100;      % [hPa]


%% ================================================================
%  CLIMATOLOGIAS DE VELOCIDADE E DIREÇÃO DO VENTO (10 m)
%  ================================================================

% Pré-alocação
clim_wspd = nan(size(clim_u10));  % velocidade [m/s]
clim_wdir = nan(size(clim_u10));  % direção [° meteorológicos]

for m = 1:nMes
    U = clim_u10(:,:,m);
    V = clim_v10(:,:,m);

    % Velocidade
    clim_wspd(:,:,m) = sqrt(U.^2 + V.^2);

    % Direção (em graus meteorológicos)
    % atan2(Y,X) → radianos; converte para graus; ajusta para 0–360° (0 = norte)
    dir = atan2(-U, -V) * (180/pi);
    dir(dir < 0) = dir(dir < 0) + 360;
    clim_wdir(:,:,m) = dir;
end



%% TESTES

% figure
% pcolor(lon,lat,t2m_C(:,:,3)')
% shading flat, shading interp
% colorbar
% colormap turbo
% ame_sul_br_pMAT([0.7 0.7 0.7],[0 0 0])
% meso_SC_pMAT([1 1 1])


%%

meses = ["Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez"];
meses_nu = ["01","02","03","04","05","06","07","08","09","10","11","12"]; 

addpath('C:\Users\augus\Downloads\reg_sul')

close all
for i = 1:12
    figure(i)
    pcolor(lon, lat, clim_swh(:,:,i)')
    shading flat; shading interp
    hold on

    %l = streamslice(lon,lat,clim_u10(:,:,i)',clim_v10(:,:,i)'); set(l,'LineWidth',.7); set(l,'Color','r');
    % hold on
    % contour(lon,lat,sp_hPa(:,:,i)','Color','r')
    % 
    c = colorbar;
    colormap parula

    ame_sul_br_pMAT([0.7 0.7 0.7],[0 0 0])
    meso_SC_pMAT([1 1 1])

    ylabel(c, 'm', 'FontSize', 12, 'FontWeight', 'bold')
    % caxis([.1876 3.2])

    title(sprintf('Significant height of combined wind waves and swell (Climatology) - %s', meses{i}), ...
          'FontSize', 14, 'FontWeight', 'bold')

    set(gca, 'Layer', 'top', 'LineWidth', 2);
    set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    set(gcf, 'WindowState', 'maximized');
    yticks(-100:2:100)
    xticks(-100:2:100)

    % 
    % cd 'C:\Users\augus\Downloads\reg_sul\climatologia\vento'
    % print(['clim_tp_' meses_nu{i}], '-djpeg', '-r600')
    % close
end

%%



%% ==== 7. Salvar resultados ====
save('climatologias_mensais_31anos.mat', ...
     'lon','lat','clim_tp','clim_t2m','clim_u10','clim_v10','clim_sp', ...
     'anos','nAnos');

fprintf('\n Climatologias mensais (31 anos) salvas com sucesso!\n');
fprintf('Cada variável é uma matriz [lon × lat × 12].\n');

%% ==== 8. (Opcional) Plotar climatologia de precipitação ====
meses = ["Jan","Fev","Mar","Abr","Mai","Jun","Jul","Ago","Set","Out","Nov","Dez"];
for m = 1:nMes
    figure('Visible','off');
    pcolor(lon, lat, clim_tp(:,:,m)'*1000); shading flat; colorbar;
    title(sprintf('Climatologia de Precipitação - %s (mm/mês)', meses(m)));
    xlabel('Longitude'); ylabel('Latitude');
    set(gca, 'YDir', 'normal');
    saveas(gcf, sprintf('clim_tp_%02d_%s.png', m, meses(m)));
end
fprintf('  Mapas mensais de precipitação salvos em mm/mês.\n');
