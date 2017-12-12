import cPickle as pickle
from os import getenv


def query():
      r = ''
      try:
            l = pickle.load(open(getenv('SUBMIT_WORKDIR')+'/submission.pkl'))
            s = l[-1]
            r += 'ClusterID'+str(s.cluster_id) + '\n'
            statii = s.query_status()
            r += 'Job summary:\n'
            for k,v in statii.iteritems():
                  r += '\t %10s : %5i\n'%(k,len(v))
            return r
      except IOError:
            print 'No job submitted yet!'

if __name__ == '__main__':
      print query()
