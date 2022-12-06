import glob

all_files   = sorted(glob.glob('./AAT_Work/*'))
stars       = all_files[:-2]
stars       = [star[11:] for star in stars]

for star in ['69830', '32147', '22049']:
    stars.remove(star)

stars = [star + '/' for star in stars]

print(stars)


with open('stars.txt', 'w') as f:
    for star in stars:
        f.write(f"{star}\n")