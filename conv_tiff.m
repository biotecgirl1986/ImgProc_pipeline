function newfilename=conv_tiff(fileName,Lx,Ly,framenum)
    fp = fopen([fileName '.tif'], 'rb');
    for cnt = 1:framenum
        imData = fread(fp, [Ly Lx], 'uint8', 0, 'ieee-be')';
        imwrite(imData,[fileName '_interleaved.tif'],'tif','compression','none','writemode','append');
    end
    fclose(fp);
    newfilename=[fileName '_interleaved'];