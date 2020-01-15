import pylab as pl
import hiv_model as hm
import tb_model as tm


fig = pl.figure(figsize=(24,16))

# Run without coinfection
hsim1 = hm.HIV_Model(hm.default_pars, coinfection=False)
tsim1 = tm.TB_Model(tm.default_pars, coinfection=False)
hsim1.run()
tsim1.run()

ax1 = pl.subplot(2,2,1)
hsim1.plot()
ax2 = pl.subplot(2,2,2)
tsim1.plot()


# Run with coinfection
hsim2 = hm.HIV_Model(hm.default_pars, coinfection=True)
tsim2 = tm.TB_Model(tm.default_pars, coinfection=True)

for t in hsim2.tvec[:-1]:
    hsim2.step()
    tsim2.step()

ax3 = pl.subplot(2,2,3)
hsim2.plot()
ax4 = pl.subplot(2,2,4)
tsim2.plot()

pl.savefig('coinfection_model.png')

print('Done.')