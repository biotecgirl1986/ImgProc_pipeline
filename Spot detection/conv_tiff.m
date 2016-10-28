function newfilename=conv_tiff(fileName,Lx,Ly,framenum)
    fp = fopen([fileName '.tif'], 'rb');
    for cnt = 1:framenum
        imData = fread(fp, [Lx Ly], 'uint16', 0, 'ieee-be')';
        imwrite(imData,[filename '_interleaved.tif'],'tif','compression','none','writemode','append');
    end
    fclose(fp);
    newfilename=[filename '_interleaved.tif'];