gasses = { 'n2', 'o2', 'ar', 'co2', 'ne', 'he', 'kr', 'xe', 'ch4', 'h2' };
mol_w = [ 28.0134, 31.9988, 39.948, 44.00995, 20.183, 4.0026, ...
          83.80, 131.30, 16.04303, 2.01594 ];
mol_frac_v = [ 0.78084, 0.209476, 0.00934, 0.000314, 0.00001818, ...
               0.00000524, 0.00000114, 0.000000087, 0.000002, 0.0000005 ];

el = { 'n', 'o', 'ar', 'c', 'ne', 'he', 'kr', 'xe', 'h' };
el_w = [mol_w(1)/2, mol_w(2)/2, mol_w(3), mol_w(4)-mol_w(2), mol_w(5), ...
        mol_w(6), mol_w(7), mol_w(8), mol_w(10)];
el_frac_v = [ mol_frac_v(1)*2, mol_frac_v(2)*2+mol_frac_v(4)*2, ...
              mol_frac_v(3), mol_frac_v(4)+mol_frac_v(9), ...
              mol_frac_v(5:8), mol_frac_v(9)*4+mol_frac_v(10)*2 ];
el_const_w = el_frac_v.*el_w;

el_frac_w = round(el_const_w/sum(el_const_w)*1e8)/1e8;
sprintf('%.8f\n',el_frac_w)

el_frac_w4 = round(el_const_w(1:4)/sum(el_const_w(1:4))*1e6)/1e6;
sprintf('%.6f\n',el_frac_w4)

el_frac_w3 = round(el_const_w(1:3)/sum(el_const_w(1:3))*1e4)/1e4;
sprintf('%.4f\n',el_frac_w3)
