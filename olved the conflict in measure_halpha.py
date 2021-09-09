[1mdiff --cc measure_halpha.py[m
[1mindex 7fbc866,981f4ca..0000000[m
[1m--- a/measure_halpha.py[m
[1m+++ b/measure_halpha.py[m
[36m@@@ -40,7 -40,7 +40,11 @@@[m [mdef make_halpha_image(halpha, galaxy)[m
      plt.colorbar()[m
      plt.xlabel("X (pix)")[m
      plt.ylabel("Y (pix)")[m
[32m++<<<<<<< HEAD[m
[32m +    #plt.savefig(f"{galaxy}_halpha.png", dpi=250)[m
[32m++=======[m
[32m+     plt.savefig(f"{galaxy}_halpha.png", dpi=250)[m
[32m++>>>>>>> c19b3236651d98a6f4f498c7816bffb56a1084c7[m
      plt.close()[m
  [m
  if __name__ == "__main__":[m
[36m@@@ -96,9 -96,6 +100,12 @@@[m
          # apertures = [CircularAperture(positions, r=r) for r in radii][m
          apertures = [][m
          plt.subplot(1, 2, 1)[m
[32m++<<<<<<< HEAD[m
[32m +        plt.title("Fluxo H-alfa")[m
[32m +        plt.xlabel('raio (px)')[m
[32m +        plt.xlabel('raio (px)')[m
[32m++=======[m
[32m++>>>>>>> c19b3236651d98a6f4f498c7816bffb56a1084c7[m
          vmin = np.percentile(halpha, 10)[m
          vmax = np.percentile(halpha, 95)[m
          plt.imshow(halpha, vmin=vmin, vmax=vmax)[m
[36m@@@ -106,37 -103,13 +113,48 @@@[m
              aperture = CircularAperture(positions, r=r)[m
              apertures.append(aperture)[m
              aperture.plot(color='r', lw=1)[m
[32m++<<<<<<< HEAD[m
[32m +        print(r)[m
          plt.subplot(1, 2, 2)[m
          phot_table = aperture_photometry(halpha, apertures, mask=mask)[m
[32m +        print(phot_table)[m
[32m++=======[m
[32m++        plt.subplot(1, 2, 2)[m
[32m++        phot_table = aperture_photometry(halpha, apertures, mask=mask)[m
[32m++>>>>>>> c19b3236651d98a6f4f498c7816bffb56a1084c7[m
          # Lendo os valores da tabela[m
          phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in[m
                  range(30)][m
          table = Table([radii, phot], names=["sma", "halpha"])[m
          table.write("photometry_halpha.fits", overwrite=True)[m
[32m++<<<<<<< HEAD[m
[32m +        plt.title("Fotometria")[m
[32m +        plt.xlabel('raio (pixel)')[m
[32m +        plt.ylabel('Brilho Superficial instrumental')[m
[32m +        plt.plot(radii, phot, "o")[m
[32m +        plt.show()[m
[32m +        #plt.savefig('halpha_photometry.png')[m
[32m +        [m
[32m +        [m
[32m +# ########Pixel Masking##############################[m
[32m +    [m
[32m +[m
[32m +    # 2) Fazer máscara das estrelas encontradas na query baseada no código[m
[32m +    # test_mascara.py[m
[32m +    # # ) Passar máscara como argumento do aperture_photometry[m
[32m +    ### Exemplo[m
[32m +   # mask = np.zeros_like(halpha)[m
[32m +    # Performing Aperture Photometry[m
[32m +    #phot_table = aperture_photometry(halpha, apertures, mask=mask)[m
[32m +    # Lendo os valores da table[m
[32m +   # phot = [float(phot_table["aperture_sum_{}".format(i)]) for i in range(30)][m
[32m +    # table = Table([radii, phot], names=["sma", "halpha"])[m
[32m +    # table.write("photometry_halpha.fits", overwrite=True)[m
[32m +    # plt.plot(radii, phot, "o")[m
[32m +    # plt.savefig('CUBE_FOTOMETRIA_Halpha.png')[m
[31m-     # plt.show()[m
[32m++    # plt.show()[m
[32m++=======[m
[32m+         plt.plot(radii, phot, "o")[m
[32m+         plt.savefig('halpha_photometry.png')[m
[31m -        plt.show()[m
[32m++        plt.show()[m
[32m++>>>>>>> c19b3236651d98a6f4f498c7816bffb56a1084c7[m
