import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs

Item {
    id: item2
    x:10
    y:10
    width: 200
    height: 400
    visible: true
    clip: true

    property bool isFixed: false

    function updateFunctionList() {
        console.log("update")
        console.log(MainCtrl.getFunctionNames())
        comboBox.model = MainCtrl.getFunctionNames()
    }

    Component.onCompleted: function() {
        close()
        MainCtrl.meshLoaded.connect(updateFunctionList)
    }

    function close(){
        if(!isFixed)
        {
            item2.anchors.top = undefined
            item2.anchors.bottom  = undefined
            height_range.to = roundButton.height
            width_range.to = roundButton.width
            widthAnim.start()
            heightAnim.start()
        }
    }

    //Component.onCompleted: close()

    SequentialAnimation on width {
        id: widthAnim
        running: false
        //PropertyAnimation { to: 36 }
        PropertyAnimation {id: width_range; to: 50 }
    }
    SequentialAnimation on height {
        id: heightAnim
        running: false
        //PropertyAnimation { to: 36 }
        PropertyAnimation {id: height_range; to: 50 }
        onFinished: function() {

            if(item2.height != roundButton.height)
            {
                item2.anchors.top = item2.parent.top
                item2.anchors.left  = item2.parent.left
                item2.anchors.bottom  = item2.parent.bottom
                 item2.anchors.bottomMargin = 10
                item2.anchors.topMargin = 10
                item2.anchors.leftMargin = 10
            }
        }
    }

    Rectangle {
        color: "#82ffffff"
        anchors.fill: parent
        radius: 7
        border.color: "#fedcdcdc"
        border.width: 1

        MouseArea {
            id: mouseArea
            anchors.fill: parent
            anchors.topMargin: 0
            anchors.bottomMargin: 0
            anchors.leftMargin: 0
            anchors.rightMargin: 0
            focus: false
            hoverEnabled: true
            drag.threshold: 10
            onClicked: {
                if(width_range.to === 200){
                }else{
                    height_range.to = appWindow.height - 20
                    width_range.to = 200
                    widthAnim.start()
                    heightAnim.start()
                    //item2.anchors.bottom = appWindow.bottom
                }
                forceActiveFocus()
            }

            RowLayout {
                id: rowLayout
                height: 32
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.top: parent.top
                anchors.topMargin: 0
                anchors.leftMargin: 0
                anchors.rightMargin: 0

                RoundButton {
                    id: roundButton1
                    opacity: 0.402
                    text: ""
                    icon.source: "../icons/settings.svg"
                    display: AbstractButton.IconOnly
                    onClicked: function() {
                        //item2.
                        height_range.to = appWindow.height - 20
                        width_range.to = 200
                        widthAnim.start()
                        heightAnim.start()
                    }
                }

                Item {
                    id: item3

                    Layout.fillHeight: true
                    Layout.fillWidth: true
                }


                RoundButton {
                    id: roundButton
                    opacity: 0.413
                    padding: 6
                    checked: false
                    checkable: true
                    highlighted: false
                    flat: false
                    display: AbstractButton.IconOnly
                    icon.source: "../icons/keep_off.svg"
                    onToggled: if(checked) {
                                   isFixed = true
                               } else if(!checked) {
                                   isFixed = false
                               }
                }



            }

            GridLayout {
                id: gridLayout
                anchors.left: parent.left
                anchors.right: parent.right
                anchors.top: rowLayout.bottom
                anchors.bottom: parent.bottom
                anchors.topMargin: 5
                columnSpacing: 5
                columns: 1

                Switch {
                    id: checkBox1
                    opacity: 0.789
                    text: qsTr("show triangle outline")
                    display: AbstractButton.TextOnly
                    Layout.fillWidth: true
                    Layout.columnSpan: 1
                    onCheckedChanged: function() {
                        MainCtrl.toogleWireframe(checked)
                    }
                }

                Switch {
                    id: checkBox
                    opacity: 0.789
                    text: qsTr("white background")
                    rotation: 0
                    Layout.fillWidth: true
                    Layout.columnSpan: 1
                }

                Switch {
                    id: checkBox2
                    opacity: 0.789
                    text: qsTr("show bounding box")
                    display: AbstractButton.TextOnly
                    Layout.columnSpan: 1
                    Layout.fillWidth: true
                }

                ColumnLayout {
                    Layout.bottomMargin: 0
                    Layout.margins: 5
                    Layout.fillWidth: true

                    Label {
                        id: label
                        opacity: 0.789
                        text: qsTr("Gridfunction")
                        Layout.bottomMargin: 0
                        Layout.margins: 5
                        Layout.fillWidth: true
                    }

                    ComboBox {
                        id: comboBox
                        opacity: 0.789
                        Layout.fillWidth: true
                        onCurrentTextChanged: function() {
                            MainCtrl.selectFunctionByName(currentText)
                        }
                    }
                }

                Item {
                    id: item1
                    Layout.fillHeight: true
                    Layout.fillWidth: true
                }
            }
        }

    }
}

/*##^##
Designer {
    D{i:0}D{i:12;invisible:true}
}
##^##*/
