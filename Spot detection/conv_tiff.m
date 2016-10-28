function newfilename=conv_tiff(fileName,framenum)
    info=imfinfo([fileName '.tif']);
    fp = fopen([fileName '.tif'], 'rb');
    fseek(fp, info.StripOffsets, 'bof');
    for cnt = 1:framenum
        imData = fread(fp, [info.Width info.Height], 'uint8', 0, 'ieee-be')';
        imwrite(imData,[filename '_interleaved.tif'],'tif','compression','none','writemode','append');
    end
    fclose(fp);
    newfilename=[fileName '_interleaved'];
    
    
    @@@ Changes to the code 
    function newfilename=conv_tiff(fileName,framenum)
    info=imfinfo([fileName '.tif']);
    fp = fopen([fileName '.tif'], 'rb');
    fseek(fp, info.StripOffsets, 'bof');
    for cnt = 1:framenum
        imData = fread(fp, [info.Width info.Height], 'uint8', 0, 'ieee-be')';
        imwrite(imData,[fileName '_interleaved.tif'],'tif','compression','none','writemode','append');
    end
    fclose(fp);
newfilename=[fileName '_interleaved'];
