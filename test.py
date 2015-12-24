from spindex4 import SpIndex4p, SpIndex4e


def main():

#     print '\n\nTESTING QUAD for points'
#     test_spindex4 = SpIndex4p((0,0,100,100))
#     print 'test_spindex4 is '
#     for cell in test_spindex4._index:
#         print cell
# #    points = {'point1': (12,15,12,15),
# #              'point2': (65,30,65,30),
# #              'point3': (40.1,85.5,40.1,85.5),
# #              'point4': (74,60,74,60),
# #              'point1': (12,15,12,15)}  # one in each original quadrant
#
#     points = {'point_a': (12,15),
#               'point_b': (49,40),
#               'point_c': (40.1,85.5),
#               'point_d': (85,60),
#               'point_e': (52,42),
#               'point_f': (12,15),
#               'point_f': (12,15)}
#
#     for point in points:
#         print 'attempting to insert point ', points[point]
#         test_spindex4.insert(points[point],point)
#         # print 'index is now'
#         # print test_spindex4._index
#
#     print '\nintersected points are:'
#     int_env = (10.0, 10.0, 70.0, 70.0)
#     int_pts = test_spindex4.intersect(int_env)
#     for pt in int_pts:
#         print pt

    # print '\n\n level 0 cells'
    # for cell in test_spindex4._index:
    #     print cell

    # print '\n\n'
    # for item in test_spindex4.walk():
    #     print '*' * item[1], item[0]

    # print '\n\nThe closest points to 10,10 within 10 are:'
    # test_nearest = test_spindex4.nearest((10.0,10.0),10.0)
    # if test_nearest:
    #     for item in test_nearest:
    #         print item

    # print '\n\n'
    # for item in set(test_spindex4.intersect((5,5,25,25))):
    #     print 'intersected item ', item
    # test_spindex4.visualize()
    #
    #
    # print '\n\n'
    # test_spindex4.print_rectangles()
    #
    # print '\n\n'
    # print test_spindex4._index
    #
    #
    # test 1 for QUAD envelopes
    print '\n\n testing quad for envelopes'
    test_spindex4e = SpIndex4e((25,25,75,75))
    # lines = {'ln1':((1,1),(3,4)),
    #          'ln2':((80,2),(4,10)),
    #          'ln3':((20,10),(10,90),(90,8),(90,90)),
    #          'ln4':((30,30),(40,40))}

    lines = {'ln1':((11.48,11.83),(30.72,12.18)),
            'ln2':((11.47,11.83),(11.48,54.71)),
            'ln3':((11.48,54.71),(60.39,54.71)),
            'ln4':((30.72,12.18),(30.72,38.83)),
            'ln5':((60.22,11.12),(30.72,12.18)),
            'ln6':((30.72,38.83),(60.22,38.83)),
            'ln7':((89.35,11.12),(60.22,11.12)),
            'ln8':((60.22,38.83),(60.22,11.12)),
            'ln9':((89.53,38.83),(89.35,11.12)),
            'ln10':((60.22,38.83),(89.53,38.83)),
            'ln11':((60.39,54.71),(60.22,38.83)),
            'ln12':((69.37,54.83),(89.53,38.83)),
            'ln13':((69.37,54.83),(60.39,54.83)),
            'ln14':((60.39,54.71),(60.46,60.85)),
            'ln15':((69.37,60.85),(69.37,54.83)),
            'ln16':((69.37,60.85),(78.83,60.85)),
            'ln17':((78.83,60.85),(78.48,67.85)),
            'ln18':((78.48,67.85),(69.41,67.85)),
            'ln19':((60.46,60.85),(69.37,60.85)),
            'ln20':((69.41,67.91),(69.37,60.85)),
            'ln21':((60.46,60.85),(60.46,67.91)),
            'ln23':((60.46,67.91),(69.41,67.91)),
            'ln24':((60.46,67.91),(60.46,74.97)),
            'ln25':((69.36,75.31),(69.41,67.91)),
            'ln26':((60.46,74.97),(69.36,75.31))}

    for line in lines:
       minx = 99999999
       miny = 99999999
       maxx = -99999999
       maxy = -99999999
       for vertex in lines[line]:
           minx = min(minx,vertex[0])
           miny = min(miny,vertex[1])
           maxx = max(maxx,vertex[0])
           maxy = max(maxy,vertex[1])
       env = (minx,miny,maxx,maxy)
       test_spindex4e.insert(env,line)

    print '\n\n'
    for cell in test_spindex4e._index:
       print cell

    print '\n\n'
    for item in test_spindex4e.walk():
       print '*' * item[1], item[0]

    print '\n\n'
    #    for item in set(test_spindex4e.intersect((20,30,40,42))):
    #    for item in set(test_spindex4e.intersect((55,70,75,94))):
    #    for item in set(test_spindex4e.intersect((10,60,30,80))):
    #    for item in set(test_spindex4e.intersect((45,35,70,60))):
    for item in set(test_spindex4e.intersect((45,35,70,60))):
       print 'intersected item ', item
    test_spindex4e.visualize()
    #
    #
    #  TEST 2 FOR QUAD
    print '\n\nTESTING QUAD'
    minx = 0
    miny = 0
    maxx = 1000
    maxy = 1000
    test4_2_bounds = (minx,miny,maxx,maxy)
    tolerance = 5.0
    test_spindex4_2 = SpIndex4p(test4_2_bounds)
    test4_2_number_of_points = 100
    test4_2_point_dict = {}
    from random import uniform
    for i in range(test4_2_number_of_points):
       x = round(uniform(minx+1,maxx-1),2)
       y = round(uniform(miny+1,maxy-1),2)
       point = (x,y)
       print 'p'+str(i), x, y
       test4_2_point_dict['p'+str(i)] = point
    for key in test4_2_point_dict:
       test_spindex4_2.insert(test4_2_point_dict[key], key)
    #
    # print '\n\n'
    # test_spindex4_2.print_rectangles()
    #
    print '\n\n'
    intersection_env = (250,250,500,500)
    print ' items intersected by ', intersection_env
    for item in set(test_spindex4_2.intersect(intersection_env)):
       print 'intersected item ', item

    test_spindex4_2.visualize()
    #
    #
    # # test on a series of random points
    # from random import uniform
    # minx = 100
    # miny = 100
    # maxx = 100000
    # maxy = 100000
    # tolerance = 5.0
    # test2_bounds = (minx,miny,maxx,maxy)
    # test2_number_of_points = 700
    # test2_spindex = SpIndex4p(test2_bounds)
    # test2_point_dict = {}
    # for i in range(test2_number_of_points):
    #    x = uniform(minx,maxx)
    #    y = uniform(miny,maxy)
    # ##        print x, y
    #    test2_point_dict['point'+str(i)] = (x,y)
    # for key in test2_point_dict:
    # ##        print key, ' is ', test2_point_dict[key]
    #    test2_spindex.insert(test2_point_dict[key], key)
    # test2_spindex.visualize()
    # test2_env_minx = uniform(minx,maxx)
    # test2_env_miny = uniform(miny,maxy)
    # test2_env_maxx = uniform(test2_env_minx,maxx)
    # test2_env_maxy = uniform(test2_env_miny,maxy)
    # test2_envelope = (test2_env_minx,
    #                  test2_env_miny,
    #                  test2_env_maxx,
    #                  test2_env_maxy)
    # print '\n the test envelope is '
    # print test2_env_minx
    # print test2_env_miny
    # print test2_env_maxx
    # print test2_env_maxy
    #
    # print '\n Items intersected by envelope ', test2_envelope, ' are:'
    # print test2_spindex.intersect(test2_envelope)
    #
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

##    test2_large_envelope = (1000,1000,90000,90000)
##    print '\n Items intersected by envelope ', test2_large_envelope, ' are:'
##    print test_spindex.intersect(test2_large_envelope)

if __name__ == '__main__':
    main()
