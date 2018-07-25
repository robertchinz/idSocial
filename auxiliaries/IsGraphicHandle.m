function tf = IsGraphicHandle(h)
        % Return TRUE if an object h is New Graphics Object
        if verLessThan('matlab','8.4')
            tf = ~isempty(h) && any(ishghandle(h));
        else
            try
        tf = ~isempty(h) && (strncmp(class(h),'matlab.ui',9) || ...
            strncmp(class(h),'matlab.graphics',15)) && ...
             any(ishghandle(h));
            catch
                keyboard
            end
        end
end