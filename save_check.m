function [data] = save_check(storage, new_data)
if isempty(storage)
    data = new_data;
else
    data = cat(3,storage,new_data);
end
end