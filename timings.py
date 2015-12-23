from spindex4 import SpIndex4p, SpIndex4e

from random import uniform
import time

def main():
    # test on a series of random points

    minx = 100
    miny = 100
    maxx = 100000
    maxy = 100000
    tolerance = 5.0
    test2_bounds = (minx,miny,maxx,maxy)
    test2_number_of_points = 100000
    test2_spindex = SpIndex4p(test2_bounds)
    test2_point_dict = {}
    for i in range(test2_number_of_points):
       x = uniform(minx,maxx)
       y = uniform(miny,maxy)
    ##        print x, y
       test2_point_dict['point'+str(i)] = (x,y)
    for key in test2_point_dict:
    ##        print key, ' is ', test2_point_dict[key]
       test2_spindex.insert(test2_point_dict[key], key)
    # test2_spindex.visualize()
    test2_env_minx = uniform(minx,maxx)
    test2_env_miny = uniform(miny,maxy)
    test2_env_maxx = uniform(test2_env_minx,maxx)
    test2_env_maxy = uniform(test2_env_miny,maxy)
    test2_envelope = (test2_env_minx,
                     test2_env_miny,
                     test2_env_maxx,
                     test2_env_maxy)
    print '\n the test envelope is '
    print test2_env_minx
    print test2_env_miny
    print test2_env_maxx
    print test2_env_maxy

    print '\n Items intersected by envelope ', test2_envelope, ' are:'
    print test2_spindex.intersect(test2_envelope)

    # for item_key in test2_point_dict:
    #    near_envelope = (test2_point_dict[item_key][0] - tolerance,
    #                     test2_point_dict[item_key][1] - tolerance,
    #                     test2_point_dict[item_key][0] + tolerance,
    #                     test2_point_dict[item_key][1] + tolerance,)
    #    intersected_items = set(test2_spindex.intersect(near_envelope))
    #    intersected_items.discard(item_key)
    #    if intersected_items:
    #       print '  features within ', tolerance, ' units of point ', item_key, ' are:'
    #       for intersectedID in intersected_items:
    #           print intersectedID




if __name__ == '__main__':

    start_time = time.time()
    start_time_display = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

    # import cProfile
    # cProfile.run('main()', sort='time')

    main()

    print
    print "*" * 80
    stop_time = time.time()
    stop_time_display = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    elapsed_time = stop_time - start_time
    print "Start Time: ", start_time_display
    print "Stop Time: ", stop_time_display
    print "Elapsed time = ",round(elapsed_time,2)," seconds"
    print "*" * 80



    # 100000 random points
    # pure python                       5.66 seconds
    # pure python compiled with Cython  3.69 seconds
    # python with some types declared   3.59 seconds
