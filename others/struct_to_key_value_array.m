function key_value_array = struct_to_key_value_array(struct_in)
    key_value_array = {}; % Initialize as empty cell array
    fields = fieldnames(struct_in);
    
    for i = 1:length(fields)
        key = fields{i};
        value = struct_in.(key);
        if isstruct(value)
            % Recursively process nested structure
            nested_key_value = struct_to_key_value_array(value);
            key_value_array = [key_value_array; nested_key_value]; % Append results
        else
            % Store key and value
            key_value_array = [key_value_array; {key, value}];
        end
    end
end