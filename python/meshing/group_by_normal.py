

def groupByNormal(normals, tollerance):
    # Erstellen Sie eine Gruppenzuordnung für Dreiecke basierend auf Normalen
    triangle_groups = {}
    triangle_groups_triangles = {}
    triangle_groupping = []
    for i, normal in enumerate(normals):
        group_id = None

        # Überprüfen Sie, ob die Normale zu einer vorhandenen Gruppe passt
        for group, group_normals in triangle_groups.items():
            for group_normal in group_normals:
                similarity = np.abs(np.dot(normal, group_normal))
                if similarity >= tolerance:
                    group_id = group
                    break
            if group_id is not None:
                break

        # Wenn keine passende Gruppe gefunden wurde, erstellen Sie eine neue Gruppe
        if group_id is None:
            group_id = len(triangle_groups)
            print("new id    : " + str(group_id))
            triangle_groups[group_id] = [normal]
            triangle_groups_triangles[group_id] = [i]
        else:
            print("add to id : " + str(group_id))
            # Fügen Sie die Normale zur entsprechenden Gruppe hinzu
            triangle_groups[group_id].append(normal)
            triangle_groups_triangles[group_id].append(i)
        
        triangle_groupping.append(group_id)