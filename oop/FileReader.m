classdef FileReader < handle
    %FILEREADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v0
        folder
        lengthscale
        coordinate
        t_index_i
        t_index_j
        vel
    end
    
    methods
        function obj = setFolder(obj, basefolder, t_index_i, t_index_j)
            obj.folder = basefolder + int2str(t_index_i) + "_" + int2str(t_index_j) + "/";
            obj.t_index_i = t_index_i;
            obj.t_index_j = t_index_j;
        end
        
        function readV0(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            v0_file = obj.folder + "v0.txt";
            try
                obj.v0 = csvread(v0_file);
            catch
                disp('no v0');
            end
        end
        
        function readInitial(obj)
            extend = "_jammed_" + int2str(obj.t_index_i) +".txt";
            coordinate_file = obj.folder + "jam" + extend;
            length_file = obj.folder + "length" + extend;
            cal_A_file = obj.folder + "calA" + extend;
            contact_file = obj.folder + "contact" + extend;

            packing_contact = dlmread(contact_file);
            packing_contact = packing_contact(end,:);
            obj.coordinate = csvread(coordinate_file);
            %length_file = folder + "length" + extend;
            obj.lengthscale = csvread(length_file);
            cal_A = csvread(cal_A_file);
        end
        
        function readMDdata(obj)
            extend1 = "_" + int2str(obj.t_index_i) + int2str(obj.t_index_j) +".txt";
            coordinate_file = obj.folder + "jam" + extend1;
            obj.coordinate = csvread(coordinate_file);
%             length_file = obj.folder + "length" + extend;
%             obj.lengthscale = csvread(length_file);
            try
                cal_A_file = obj.folder + "calA" + extend1;
                cal_A = csvread(cal_A_file);
            catch
                disp('no calA file')
            end
            try
                contact_file = obj.folder + "contact" + extend1;
                contact = dlmread(contact_file);
            catch
                disp('no contact file')
            end
            v_file = obj.folder + "v" + extend1;
            obj.vel = csvread(v_file);
        end
    end
end

