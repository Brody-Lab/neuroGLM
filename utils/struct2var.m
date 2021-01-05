function varargout = struct2var(struct,fields,cat_dim,add_string)
    if nargin<4
        add_string='';       
        if nargin<3
            cat_dim=1;
            if nargin<2 
                fields='all';
            else
                if ischar(fields)
                    fields={fields};
                end
            end
        end
    end
    if strcmp(fields,'all')
        fields=fieldnames(struct);
    end
    if ~isstruct(struct) 
        error('First input must be scalar structure with fields');
    end
    if ~numel(fields)
        error('No fields');
    end
    if nargout==0
        for i=1:length(fields)
            if ischar(struct(1).(fields{i}))
                assignin('caller',[fields{i},add_string],{struct.(fields{i})});                
            else
                assignin('caller',[fields{i},add_string],cat_nans(cat_dim,struct.(fields{i})));
            end
        end
    else
        if nargout~=numel(fields)
            error('Number of outputs must match number of fields.');
        end
        for i=nargout:-1:1
            varargout{i}=cat_nans(cat_dim,struct.(fields{i}));
        end
    end
end

function d= cat_nans(dim,varargin)
    nans = cellfun(@isempty,varargin);
    varargin(nans)={NaN};
    d=cat(dim,varargin{:});

end