import pandas as pd



a=pd.DataFrame(data=[[1,2,5,8],[3,7,9,10],[4,10,61,7]],index=['a','b','c'],columns=['cc','bb','bbe','dde'])
# print(a)
c = a.loc[['a','b']].sum(axis=0)
c=pd.Series.to_frame(c).T
c.index=['rr']
# print(ff)
e=a.drop(['a','b'])
# print(e)
print(c)
d=pd.concat([e,c])
print(d)



a='ABC CIN'
if 'ABC' not in a :
    print('yes')