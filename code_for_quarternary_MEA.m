%MATLAB code for generating quaternary medium entropy alloys 

close all;
clear all
clc; 

R = 8.3144598;

% Reading data from tables and converting to array
delH_mix = table2array(readtable('Mixing enthalpy.xlsx','ReadRowNames',true));
Element_prop = table2array(readtable('Element properties.xlsx','ReadRowNames',true)); %'ReadVariableNames',false

%Elements = table2array(Element_prop_table(:,1));
Density_element = Element_prop(:,1);
Atomic_radius_element = Element_prop(:,2);
Elastic_modulus_element = Element_prop(:,3);
Cost_element = Element_prop(:,4);
Melting_temp_element = Element_prop(:,5);
Mol_weight = Element_prop(:,6);

% Array of 20 bio-metals
elements_ar = {'Mg','Zn','Ti','Zr','Fe','Ca','Cu','Sr','Nb','Hf','Ta','Re','Cr','Mo','Co','Y','Ag','Bi','Mn','Yb'};
Len_elements_ar = length(elements_ar);
element_conc = [0.25 0.25 0.25 0.25]; %Equiatomic composition
%s = size(element_conc);

% Variable initialization
N_alloys = nchoosek(Len_elements_ar,4);  %combination (20C4=4845), 
k=0;
Alloy = cell(N_alloys,1);
Mixing_enthalpy = zeros(N_alloys,1);
Tm = zeros(N_alloys,1);
Omega = zeros(N_alloys,1);
Delta = zeros(N_alloys,1);
Cost = zeros(N_alloys,1);
Density = zeros(N_alloys,1);
Youngs_modulus = zeros(N_alloys,1);


for a = 1:Len_elements_ar
    element_a = elements_ar(a);
    for b = a+1:Len_elements_ar
        element_b = elements_ar(b);
        for c = b+1:Len_elements_ar
            element_c = elements_ar(c);
            for d = c+1:Len_elements_ar
                element_d = elements_ar(d);
                
                
                e = cell2mat({find(strcmp(elements_ar,element_a)),find(strcmp(elements_ar,element_b)),find(strcmp(elements_ar,element_c)),find(strcmp(elements_ar,element_d))});
                n = length(e);
                alloy = {strjoin([element_a,element_b,element_c,element_d], '')};
                
                
                density_avg = 0;
                radius_avg = 0;
                Tm_avg = 0;
                modulus_avg = 0;
                cost_avg = 0;
                delta_Hmix = 0;
                delta = 0;
                weight = 0;
                del_Smix = 0;
                d_numerator = 0;
                d_denominator = 0;
                density = zeros(n,1);
                radius = zeros(n,1);
                tm = zeros(n,1);
                cost = zeros(n,1);
                modulus = zeros(n,1);
                m_weight = zeros(n,1);
                
                for i = 1:n
                    density(i)= Density_element(e(i));
                    radius(i)= Atomic_radius_element(e(i));
                    tm(i)= Melting_temp_element(e(i));
                    cost(i) = Cost_element(e(i));
                    modulus(i) = Elastic_modulus_element(e(i));
                    m_weight(i) = Mol_weight(e(i));
                    
                    radius_avg = radius(i)*element_conc(i)+radius_avg;
                    Tm_avg = tm(i)*element_conc(i)+Tm_avg;
                    modulus_avg = modulus(i)*element_conc(i)+modulus_avg;
                    weight = m_weight(i)*element_conc(i)+weight;
                    d_numerator = element_conc(i)*m_weight(i)+d_numerator;
                    d_denominator = (element_conc(i)*m_weight(i))/density(i)+d_denominator;
                end
                density_avg = d_numerator/d_denominator;

                %Calculation of mixing enthalpy
                for i = 1:n
                    for j = i+1:n
                        delta_Hmix = 4*delH_mix(e(i),e(j))*element_conc(i)*element_conc(j)+delta_Hmix;
                    end
                end
                
                %Calculation of atomic size mismatch (delta)
                for i = 1:n
                    delta = element_conc(i)*(1-(radius(i)/radius_avg))^2+delta;
                    cost_avg = cost(i)*element_conc(i)*m_weight(i)/weight+cost_avg;
                end
                delta = sqrt(delta)*100;
                
                %Calculation of mixing entropy
                for i = 1:n
                    del_Smix = -R*element_conc(i)*log(element_conc(i))+del_Smix;
                end
                
                k=k+1;
                Mixing_enthalpy(k,1) = delta_Hmix;
                Tm(k,1) = Tm_avg;
                Density(k,1) = density_avg;
                Omega(k,1) = (Tm_avg*del_Smix)/abs(delta_Hmix*1000);
                Delta(k,1) = delta;
                Youngs_modulus(k,1) =  modulus_avg;
                Cost(k,1) = cost_avg;
                Alloy(k,1) = alloy;
                        
            end
        end
    end
end

%Keeping important variables
clearvars -except Alloy Mixing_enthalpy Tm Omega Delta Cost Density Youngs_modulus

%Creating Table containing all possible quarternary alloys
T = table(Alloy,Tm,Density,Cost,Youngs_modulus,Mixing_enthalpy,Omega,Delta);
display(T);

%Creating an excel file containing the alloy information
filename = 'Quarternary_alloys.xlsx';
writetable(T,filename,'Sheet',1,'Range','A1');
