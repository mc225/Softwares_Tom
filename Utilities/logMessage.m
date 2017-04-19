% message=logMessage(message,values,recipientEmailToAddresses,attachments)
%
% Logs a message to the command line and to a file called 'log.txt' in the
% current folder with a timestamp prefix. Optionally, it can be sent by e-mail
% as well with (optional) attachments.
% The message is formatted using escape characters when a second argument
% is provided,  containing either the substitions, or an empty list [].
%
% Input arguments:
%    message: a text string, which may be formatted using the sprintf
%             syntax (see help sprintf)
%    values: optional values that can be used by a formatted message.
%    recipientEmailToAddresses: optional email address notify, or a cell array thereof.
%    attachments: optional string or cell array of strings with filenames
%    to attachments to be included with the email
%
% Output arguments:
%    message: the formatted message
%
% Examples:
%    logMessage('This is a simple message.');
%    2013-01-26 17:08:56.589| This is a simple message
%
%    logMessage('This is a formatted message:\n pi=%f, e=%f!',[pi exp(1)]);
%    2013-01-26 17:14:9.627| This is a formatted message:
%     pi=3.141593, e=2.718282!
%
%    logMessage('This is a formatted message\nwithout substitutions.',[]);
%    2013-01-26 17:14:56.122| This is a formatted message
%    without substitutions.
%
%    logMessage('This message is e-mailed to me@here.com',[],'me@here.com');
%    2013-01-26 17:16:46.277| This message is e-mailed to me@here.com
%
%    logMessage('And this one e-mails an attachment with it.',[],'me@here.com','output.mat');
%    2013-01-26 17:18:46.409| And this one e-mails an attachment with it.
%
function message=logMessage(message,values,recipientEmailToAddresses,attachments)
    timeNow=clock();
    if (nargin<1 || isempty(message))
        message='';
    end
    if (isnumeric(message))
        message=num2str(message);
    end
    if (nargin>1),
        message=sprintf(message,values);
    end
    if (nargin<3)
        recipientEmailToAddresses=[];
    end
    if (nargin<4)
        attachments={};
    else
        if (ischar(attachments))
            attachments={attachments};
        end
    end
    %Prepare the message
    synopsis=message;
    message=sprintf('%04.0f-%02.0f-%02.0f %02.0f:%02.0f:%06.3f| %s',timeNow(1),timeNow(2),timeNow(3),timeNow(4),timeNow(5),timeNow(6),message);
    disp(message);
    %Write the message to file
    fid = fopen('log.txt','a');
    if (fid>0)
        fprintf(fid,'%s\n',message);
        fclose(fid);
    end
    %Also e-mail this message if it was requested
    if (~isempty(recipientEmailToAddresses) && (~iscell(recipientEmailToAddresses) || ~isempty(recipientEmailToAddresses{1})))
        if (~ispref('Internet','SMTP_Server') || isempty(getpref('Internet','SMTP_Server')))
            setpref('Internet','SMTP_Server','mailhost.st-andrews.ac.uk');
        end
        if (~ispref('Internet','E_mail') || isempty(getpref('Internet','E_mail')))
            setpref('Internet','E_mail','tv2@st-andrews.ac.uk');
        end
        synopsis=regexprep(synopsis,'\s+',' ');
        if (length(synopsis)>200)
            synopsis=[synopsis(1:150) '...' synopsis(end-46:end)];
        end
        try
            existingAttachments={};
            for attachmentIdx=1:length(attachments)
                attachment=attachments{attachmentIdx};
                if (ischar(attachment))
                    if (exist(attachment,'file'))
                        existingAttachments{end+1}=attachment;
                    else
                        message=[message sprintf('\n COULD NOT ATTACH FILE %s!',attachment)];
                    end
                else
                    message=[message sprintf('\n SPECIFIED ELEMENT IS NOT A FILENAME!')];
                end
            end
            try
                sendmail(recipientEmailToAddresses,synopsis,message,existingAttachments);
            catch Exc
                % Try without attachements
                sendmail(recipientEmailToAddresses,synopsis,[message sprintf('\n\n\n SENDING WITH ATTACHMENT FAILED!')]);
            end
        catch Exc
            Exc
            logMessage('Could not send as e-mail: %s',Exc.message);
        end
    end
end