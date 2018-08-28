
exhaust_filename = 'Exhaust_Lift_Profile.csv';
exhaust_data = csvread(exhaust_filename, 1, 0);

exhaust_opens = exhaust_data(1, 3);
exhaust_closes = exhaust_data(1, 4);
exhaust_lift_profile = exhaust_data(:, 1:2);


inlet_filename = 'Inlet_Lift_Profile.csv';
inlet_data = csvread(inlet_filename, 1, 0);

inlet_opens = inlet_data(1, 3);
inlet_closes = inlet_data(1, 4);
inlet_lift_profile = inlet_data(:, 1:2);